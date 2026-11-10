module Placement

    using Logging
    using Distributed
    using Distances

    using ..FileHandling: binary_path,group_meta_path,load_groups_as_dictionary,write_group_results,mv_milk_files,metadata_path,load_milk_binaries
    using ..IncrementalFileHandling: final_iteration_path,write_final_iteration
    using ..Linkage: squared_ward_cost
    using ..PairwiseComparisons: compute_pairwise_distance_matrix,compute_euclidean_distance
    using ..GroupStratification: ward_stratification

    export BackboneIteration,
           load_backbone,
           assign_to_backbone,
           merge_leftovers,
           place_stage,
           place_stages

    """
    One iteration of the backbone (the finals). Node j is column j of centroids.
    members[j]   indices of node j's members in the iteration below; empty at iteration 0 (its members are cells)
    threshold    Ward distance threshold t of the iteration (two objects merge when their Ward cost <= t^2)
    backbone[k+1] holds iteration k
    """
    struct BackboneIteration
        ids::Vector{String}
        centroids::Matrix{Float32}
        sizes::Vector{Int}
        spreads::Vector{Float64}
        radii::Vector{Float64}
        stages::Vector{Int}
        members::Vector{Vector{Int}}
        threshold::Float64
    end

    function load_backbone_iteration(ids_path,previous)
        ids = readlines(ids_path)
        n = length(ids)
        d = div(filesize(binary_path(ids_path)),4*n) # Float32 per centroid
        centroids = Matrix{Float32}(undef,d,n)
        open(binary_path(ids_path)) do io
            read!(io,centroids)
        end

        sizes = zeros(Int,n)
        spreads = zeros(Float64,n)
        radii = zeros(Float64,n)
        stages = zeros(Int,n)
        for (j,line) in enumerate(eachline(group_meta_path(ids_path)))
            fields = split(line,'\t') # id, size, spread, radius, stage (same order as .ids)
            sizes[j] = parse(Int,fields[2])
            spreads[j] = parse(Float64,fields[3])
            radii[j] = parse(Float64,fields[4])
            stages[j] = parse(Int,fields[5])
        end

        groups_dict = load_groups_as_dictionary(replace(ids_path,r"\.ids$" => ".groups.jsonl.gz"))
        threshold = parse(Float64,groups_dict["threshold"])
        members = Vector{Vector{Int}}(undef,n)
        if isnothing(previous)
            for j in 1:n
                members[j] = Int[]
            end
        else
            previous_index = Dict{String,Int}()
            for (j,id) in enumerate(previous.ids)
                previous_index[id] = j
            end
            for j in 1:n
                members[j] = Int[]
                for member in groups_dict["groups"][ids[j]]
                    push!(members[j],previous_index[member])
                end
            end
        end
        return BackboneIteration(ids,centroids,sizes,spreads,radii,stages,members,threshold)
    end

    function load_backbone(output_dir,label)
        """ all final iterations, in order: backbone[k+1] = iteration k """
        backbone = BackboneIteration[]
        iteration = 0
        while isfile(final_iteration_path(output_dir,label,iteration))
            previous = iteration == 0 ? nothing : backbone[end]
            push!(backbone,load_backbone_iteration(final_iteration_path(output_dir,label,iteration),previous))
            iteration += 1
        end
        if isempty(backbone)
            error("No backbone (final) iterations for $(label) in $(output_dir)")
        end
        @info "Loaded backbone: $(length(backbone)) iterations; nodes per iteration: $([length(it.ids) for it in backbone])"
        return backbone
    end

    function _assign_chunk(node_centroids,node_sizes,X_query,query_sizes)
        """ one chunk of query objects: each one's cheapest node and its Ward cost (runs on a worker) """
        dist_squared = pairwise(SqEuclidean(),node_centroids,Float64.(X_query);dims=2) # nodes x queries
        n_nodes,n_queries = size(dist_squared)
        best_node = zeros(Int,n_queries)
        best_cost = fill(Inf,n_queries)
        for i in 1:n_queries
            for j in 1:n_nodes
                cost = squared_ward_cost(dist_squared[j,i],query_sizes[i],node_sizes[j])
                if cost < best_cost[i]
                    best_cost[i] = cost
                    best_node[i] = j
                end
            end
        end
        return best_node,best_cost
    end

    function assign_to_backbone(iteration,X_query,query_sizes;chunk_size=1000)
        """
        Each query object (column i of X_query, size query_sizes[i]) against every node of one backbone iteration.
        Returns, per query object, the cheapest node by Ward cost (0 if that cost exceeds t^2: a leftover) and the cost.
        Chunks of chunk_size query objects are spread over workers; each holds a nodes x chunk_size Float64 block
        (e.g., 100k nodes x 1000 queries = 800 MB).
        """
        node_centroids = Float64.(iteration.centroids)
        n_queries = size(X_query,2)
        chunks = UnitRange{Int}[]
        for start in 1:chunk_size:n_queries
            push!(chunks,start:min(start+chunk_size-1,n_queries))
        end
        inputs = [(X_query[:,chunk],query_sizes[chunk]) for chunk in chunks]
        pool = CachingPool(workers()) # node centroids are sent to each worker once, not per chunk
        results = pmap(input -> _assign_chunk(node_centroids,iteration.sizes,input[1],input[2]),pool,inputs)

        best_node = zeros(Int,n_queries)
        best_cost = zeros(Float64,n_queries)
        for (chunk,(chunk_node,chunk_cost)) in zip(chunks,results)
            best_node[chunk] = chunk_node
            best_cost[chunk] = chunk_cost
        end

        max_cost = iteration.threshold^2*(1+1e-5) # same tolerance as ward_groups (linkage.jl)
        for i in 1:n_queries
            if best_cost[i] > max_cost
                best_node[i] = 0 # leftover
            end
        end
        return best_node,best_cost
    end


    function merge_leftovers(leftover_dict,leftover_size_dict,threshold)
        """
        Leftover query objects at one level (id => centroid, id => size), merged among themselves by Ward at the
        backbone's frozen threshold t, with MILK's own ward_stratification (same costs, same nearest-neighbour-chain merging).
        Returns the new nodes: id => centroid, id => size, id => members (ids of the leftovers it merged).
        A leftover with no partner within t becomes a node of its own.
        Memory: an m x m distance matrix for m leftovers.
        """
        if isempty(leftover_dict)
            return Dict{String,Vector{Float32}}(),Dict{String,Int}(),Dict{String,Vector{String}}()
        end
        D,index_map,_,_ = compute_pairwise_distance_matrix(leftover_dict,nothing,compute_euclidean_distance,false)
        node_dict,members_dict,_,_ = ward_stratification(leftover_dict,leftover_size_dict,D,index_map,threshold,nothing)

        node_size_dict = Dict{String,Int}()
        for (id,members) in members_dict
            node_size_dict[id] = 0
            for member in members
                node_size_dict[id] += leftover_size_dict[member]
            end
        end
        return node_dict,node_size_dict,members_dict
    end


    function place_stage(backbone,query_dict;chunk_size=10000)
        """
        Level-by-level placement of one stage's cells (query_dict: id => vector) onto the backbone.
        At iteration k (level k+1), each query object either joins its cheapest backbone node (Ward cost <= t_k^2)
        or is a leftover; leftovers are merged among themselves at t_k into new nodes, which are the query objects at k+1.
        New nodes still left after the top iteration attach to the root.
        Returns, per iteration k (index k+1):
            joins[k+1]       query object id => id of the backbone node it joins
            new_nodes[k+1]   (id => centroid, id => size, id => members) of the nodes formed from leftovers
        """
        query_size_dict = Dict{String,Int}()
        for id in keys(query_dict)
            query_size_dict[id] = 1 # cells
        end

        joins = Dict{String,String}[]
        new_nodes = []
        for k in 0:length(backbone)-1
            iteration = backbone[k+1]
            ids = sort(collect(keys(query_dict)))
            X_query = stack(query_dict[id] for id in ids)
            query_sizes = [query_size_dict[id] for id in ids]
            best_node,_ = assign_to_backbone(iteration,X_query,query_sizes;chunk_size=chunk_size)

            iteration_joins = Dict{String,String}()
            leftover_dict = Dict{String,Vector{Float32}}()
            leftover_size_dict = Dict{String,Int}()
            for (i,id) in enumerate(ids)
                if best_node[i] > 0
                    iteration_joins[id] = iteration.ids[best_node[i]]
                else
                    leftover_dict[id] = query_dict[id]
                    leftover_size_dict[id] = query_size_dict[id]
                end
            end
            node_dict,node_size_dict,members_dict = merge_leftovers(leftover_dict,leftover_size_dict,iteration.threshold)
            push!(joins,iteration_joins)
            push!(new_nodes,(node_dict,node_size_dict,members_dict))
            @info "\tIteration $k: $(length(ids)) query objects; $(length(iteration_joins)) join the backbone; $(length(leftover_dict)) leftovers -> $(length(node_dict)) new nodes"

            if isempty(node_dict) # everything placed
                break
            end
            query_dict = node_dict
            query_size_dict = node_size_dict
        end
        return joins,new_nodes
    end


    function placement_groups_path(stage_dir,label,iteration)
        return joinpath(stage_dir,"$(label).iteration_$(lpad(iteration,8,'0')).placement.groups.jsonl.gz")
    end

    function write_placement_groups(backbone,joins,new_nodes,output_dir,stage_dir,label)
        """
        Updated grouping at every backbone iteration k after placing one stage, written as a MILK groups file
        (metadata holds the frozen threshold t_k): stage_dir/<label>.iteration_k.placement.groups.jsonl.gz
            backbone node   its old members plus the query objects that joined it
            new node        the leftovers merged into it
        Iterations above where placement finished keep their members; sizes are recomputed from members at every iteration.
        joins and new_nodes as returned by place_stage.
        """
        previous_size_dict = Dict{String,Int}() # sizes of the iteration below's nodes
        for k in 0:length(backbone)-1
            iteration = backbone[k+1]
            old_groups_path = replace(final_iteration_path(output_dir,label,k),r"\.ids$" => ".groups.jsonl.gz")
            groups = load_groups_as_dictionary(old_groups_path)["groups"]
            if k+1 <= length(joins)
                for (query_id,node_id) in joins[k+1]
                    push!(groups[node_id],query_id)
                end
                for (id,members) in new_nodes[k+1][3]
                    groups[id] = members
                end
            end

            total_size_dict = Dict{String,Int}()
            direct_groupsize_dict = Dict{String,Int}()
            distances_dict = Dict{String,Vector{Float32}}()
            specificity_dict = Dict{String,Vector{Int32}}()
            for (id,members) in groups
                total_size_dict[id] = 0
                for member in members
                    if k == 0
                        total_size_dict[id] += 1 # cells
                    else
                        total_size_dict[id] += previous_size_dict[member]
                    end
                end
                direct_groupsize_dict[id] = length(members)
                distances_dict[id] = Float32[] # not computed in placement; recomputed exactly in the finals
                specificity_dict[id] = Int32[]
            end

            write_group_results(
                path=placement_groups_path(stage_dir,label,k),
                label="$(label).iteration_$(lpad(k,8,'0'))",
                stage="placement",
                n_cached=length(groups),
                n_input_objects=sum(values(total_size_dict)),
                n_groups=length(groups),
                groups=groups,
                direct_groupsize_dict=direct_groupsize_dict,
                total_groupsize_dict=total_size_dict,
                distances_dict=distances_dict,
                specificity_dict=specificity_dict,
                threshold=iteration.threshold,
                n_comparisons=0
            )
            previous_size_dict = total_size_dict
        end
    end


    function reconcile_backbone(backbone,output_dir,stage_dir,label,stage,cell_ids_paths)
        """
        New finals after placing one stage: every iteration rebuilt bottom-up from the placement groups
        (write_final_iteration recomputes centroid, size, spread, radius; existing nodes keep their stage, new nodes get `stage`).
        Written to <label>.iteration_k.final.tmp.* first and moved over the old finals only once every iteration is written,
        so a failure leaves the old backbone intact. The placement groups are deleted afterwards.
        cell_ids_paths: input binaries of every stage placed so far (iteration 0's members are cells)
        """
        tmp_ids_paths = String[]
        previous_ids_paths = cell_ids_paths
        for k in 0:length(backbone)-1
            iteration = backbone[k+1]
            node_stage_dict = Dict{String,Int}()
            for (j,id) in enumerate(iteration.ids)
                node_stage_dict[id] = iteration.stages[j]
            end
            tmp_ids_path = replace(final_iteration_path(output_dir,label,k),r"\.ids$" => ".tmp.ids")
            write_final_iteration(previous_ids_paths,placement_groups_path(stage_dir,label,k),tmp_ids_path,stage;node_stage_dict=node_stage_dict)
            push!(tmp_ids_paths,tmp_ids_path)
            previous_ids_paths = [tmp_ids_path]
        end

        for (k,tmp_ids_path) in zip(0:length(backbone)-1,tmp_ids_paths)
            ids_path = final_iteration_path(output_dir,label,k)
            tmp_groups_path = replace(tmp_ids_path,r"\.ids$" => ".groups.jsonl.gz")
            groups_path = replace(ids_path,r"\.ids$" => ".groups.jsonl.gz")
            mv_milk_files(tmp_ids_path,ids_path;force=true) # .ids/.bin/.group.tsv
            mv(tmp_groups_path,groups_path;force=true)
            mv(metadata_path(tmp_groups_path),metadata_path(groups_path);force=true)
            rm(placement_groups_path(stage_dir,label,k))
            rm(metadata_path(placement_groups_path(stage_dir,label,k)))
        end
        @info "Reconciled backbone: $(length(backbone)) iterations written to $(output_dir)"
    end


    function place_stages(stages_path,output_dir,label;chunk_size=10000)
        """
        Every stage after the first, in order: load the backbone (finals), place the stage's cells, write the
        placement groups to output_dir/stage_XXXXXXXX/, then reconcile them into new finals.
        stages_path: <file_label>.stages.tsv (stage, input ids path, number of cells); its first line is the stage-1 backbone
        """
        cell_ids_paths = String[]
        for line in eachline(stages_path)
            fields = split(line,'\t')
            stage = parse(Int,fields[1])
            stage_ids_path = String(fields[2])
            if isempty(cell_ids_paths)
                push!(cell_ids_paths,stage_ids_path) # stage 1: already the backbone
                continue
            end

            @info "Placing stage $(stage) ($(fields[3]) cells)"
            backbone = load_backbone(output_dir,label)
            query_dict,_ = load_milk_binaries(stage_ids_path)
            stage_dir = joinpath(output_dir,"stage_$(lpad(stage,8,'0'))")
            mkpath(stage_dir)
            joins,new_nodes = place_stage(backbone,query_dict;chunk_size=chunk_size)
            write_placement_groups(backbone,joins,new_nodes,output_dir,stage_dir,label)
            push!(cell_ids_paths,stage_ids_path)
            reconcile_backbone(backbone,output_dir,stage_dir,label,stage,cell_ids_paths)
        end
    end

end