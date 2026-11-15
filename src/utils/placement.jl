module Placement

    using Glob
    using Logging
    using Random
    using Distributed
    using Distances

    using ..FileHandling: load_groups_as_dictionary,mv_milk_files,metadata_path,load_milk_binary_matrix,
                          load_milk_group_information,get_threshold,get_iteration,jsonl_path,tree_level_path,write_tree_level
    using ..Partitioning: generalized_kmeans
    using ..Linkage: squared_ward_cost

    export ReferenceLevel,
           load_reference_level,
           get_reference_pathlist,
           greedy_placement,
           place_query_objects_to_reference,
           recompute_backbone,
           placement

    struct ReferenceLevel
        label::String
        ids::Vector{String}
        X::Matrix{Float32}
        sizes::Vector{Int}
        stages::Vector{Int}
        threshold::Float64
    end

    function load_reference_level(ids_path)
        ids,X = load_milk_binary_matrix(ids_path)
        info_dict = load_milk_group_information(ids_path)
        sizes = [info_dict.size[id] for id in ids]
        stages = [info_dict.stage[id] for id in ids]
        threshold = get_threshold(ids_path)
        return ReferenceLevel(ids_path,ids,X,sizes,stages,threshold)
    end

    function get_reference_pathlist(output_dir,label)
        groups_pathlist = sort(glob("$(label).iteration_????????.groups.jsonl.gz",output_dir))
        if isempty(groups_pathlist)
            error("No groups files for $(label) in $(output_dir)")
        end
        reference_pathlist = [replace(path,r"\.groups\.jsonl\.gz$" => ".ids") for path in groups_pathlist]
        for (index,path) in enumerate(reference_pathlist)
            if get_iteration(path) != index-1
                error("Iterations of $(label) are not consecutive: $(basename(path)) at position $(index-1)")
            end
        end
        n_root = parse(Int,split(readline(metadata_path(groups_pathlist[end])),'\t')[6]) # n_groups
        if n_root != 1
            error("Top iteration $(basename(groups_pathlist[end])) has $(n_root) nodes; expected a single root")
        end
        @info "Reference: $(length(reference_pathlist)) iterations"
        return reference_pathlist
    end

    function compute_reference_level(ids_path,previous_level)
        """ reference level computed in memory from its groups and the level below (no binaries on disk) """
        groups = load_groups_as_dictionary(jsonl_path(ids_path))["groups"]
        previous_index = Dict(id => i for (i,id) in enumerate(previous_level.ids))
        ids = collect(keys(groups))
        X = Matrix{Float32}(undef,size(previous_level.X,1),length(ids))
        sizes = zeros(Int,length(ids))
        for (j,id) in enumerate(ids)
            members = [previous_index[member] for member in groups[id]]
            sizes[j] = sum(previous_level.sizes[members])
            X[:,j] = _weighted_centroid(previous_level.X[:,members],previous_level.sizes[members])
        end
        return ReferenceLevel(ids_path,ids,X,sizes,zeros(Int,length(ids)),get_threshold(ids_path))
    end

    function _identify_min_cost_mappings(X_reference,reference_sizes,X_query,query_sizes;skip_self=false,partition_size=10000)
        """
        Each query object's cheapest reference node
        skip_self necessary when reference == query (partner search)
        Returns min_nodes (0 if there is no other node), min_costs
        """
        n_reference = size(X_reference,2)
        n_query = size(X_query,2)

        min_nodes = zeros(Int,n_query)
        min_costs = fill(Inf,n_query)
        for partition in Iterators.partition(1:n_reference,partition_size)
            D = pairwise(SqEuclidean(),view(X_reference,:,partition),X_query;dims=2)
            D .= squared_ward_cost.(D,reshape(query_sizes,1,:),reference_sizes[partition])
            if skip_self
                for (row,i) in enumerate(partition) # node i (= query i) sits at row `row` of D
                    D[row,i] = Inf
                end
            end
            for i in 1:n_query
                row = argmin(view(D,:,i))
                if D[row,i] < min_costs[i]
                    min_costs[i] = D[row,i]
                    min_nodes[i] = partition[row]
                end
            end
        end
        return min_nodes,min_costs
    end

    function _weighted_centroid(X,sizes)
        return X*sizes/sum(sizes)
    end

    function greedy_placement(reference,X_query,query_sizes)
        max_cost = reference.threshold^2*(1+1e-5)
        n_queries = size(X_query,2)
        sizes = copy(query_sizes)
        members = [[i] for i in 1:n_queries]
        remaining = trues(n_queries) # query objects that need to be placed
        query_mappings = zeros(Int,n_queries) # pos: query objects (leaf or group); value is index into reference.ids

        reference_nodes,reference_costs = _identify_min_cost_mappings(reference.X,reference.sizes,X_query,sizes)
        query_nodes,query_costs = _identify_min_cost_mappings(X_query,sizes,X_query,sizes;skip_self=true)

        function _reference_node_cost_update!(i)
            node,cost = _identify_min_cost_mappings(reference.X,reference.sizes,X_query[:,i:i],sizes[i:i])
            reference_nodes[i] = node[1]
            reference_costs[i] = cost[1]
        end

        function _query_node_cost_update!(i)
            candidates = [j for j in 1:n_queries if remaining[j] && j != i]
            node,cost = _identify_min_cost_mappings(view(X_query,:,candidates),sizes[candidates],X_query[:,i:i],sizes[i:i])
            query_nodes[i] = node[1] == 0 ? 0 : candidates[node[1]]
            query_costs[i] = cost[1]
        end

        function _remove_query!(i) # removed objects (mapped, or absorbed by a merge) carry cost Inf, so argmin never picks them
            remaining[i] = false
            query_costs[i] = Inf
            reference_costs[i] = Inf
        end

        while true
            query_merge_object = argmin(query_costs)
            reference_merge_object = argmin(reference_costs)
            query_merge_cost = query_costs[query_merge_object]
            reference_merge_cost = reference_costs[reference_merge_object]
            if min(query_merge_cost,reference_merge_cost) > max_cost
                break
            end

            # reference merge: every member joins the reference node
            if reference_merge_cost <= query_merge_cost
                i = reference_merge_object
                query_mappings[members[i]] .= reference_nodes[i]
                _remove_query!(i)
                for k in 1:n_queries
                    if remaining[k] && query_nodes[k] == i
                        _query_node_cost_update!(k)
                    end
                end
            else # query merge: j is absorbed into i; i is re-assessed against the reference and the queries
                i = query_merge_object
                j = query_nodes[i]
                X_query[:,i] = _weighted_centroid(X_query[:,[i,j]],sizes[[i,j]])
                sizes[i] += sizes[j]
                append!(members[i],members[j])
                _remove_query!(j)
                _reference_node_cost_update!(i)
                for k in 1:n_queries
                    if remaining[k] && (k == i || query_nodes[k] == i || query_nodes[k] == j)
                        _query_node_cost_update!(k)
                    end
                end
            end
        end
        unmapped_groups = [members[i] for i in 1:n_queries if remaining[i]]
        return query_mappings,unmapped_groups
    end

    function partition_queries(X_query,partition_size,rng)
        n = size(X_query,2)
        k = cld(n,partition_size)
        if k == 1
            return [collect(1:n)]
        end
        _,assignments = generalized_kmeans(X_query,k,SqEuclidean(),rng) # similarity-aware partitioning
        partition_dict = Dict{Int,Vector{Int}}() # cluster label => query objects
        for (i,assignment) in enumerate(assignments)
            push!(get!(partition_dict,assignment,Int[]),i)
        end
        return collect(values(partition_dict))
    end

    function place_query_objects_to_reference_level(reference,X_query,query_sizes,rng;partition_size=10000)
        partitions = partition_queries(X_query,partition_size,rng)
        input_list = [ (X_query[:,partition],query_sizes[partition]) for partition in partitions ]
        pool = CachingPool(workers()) # the reference level is sent to each worker once
        results = pmap(input -> greedy_placement(reference,input[1],input[2]),pool,input_list)

        n = size(X_query,2)
        query_mappings = zeros(Int,n)
        unmapped_groups = Vector{Vector{Int}}()
        for (partition,(_query_mappings,_unmapped_groups)) in zip(partitions,results)
            query_mappings[partition] = _query_mappings
            for _unmapped_group in _unmapped_groups
                push!(unmapped_groups,partition[_unmapped_group])
            end
        end
        return query_mappings,unmapped_groups
    end

    function place_query_objects_to_reference(reference_pathlist,reference_leaves_pathlist,query_ids,X_query,query_label,rng;partition_size=10000)
        """
        Bottom-up placement of all query cells. At each iteration, every query object joins a reference node (mapped) or ends up
        in an unmapped group, forming a new node holding only query objects; the unmapped groups are the query objects of the next
        iteration. Whatever is left at the top iteration joins the root.

        Returns one entry per reference iteration (index = iteration+1; empty once placement has stopped):
            query_mappings_list[index]   query object id => reference node id it joined
            unmapped_groups_list[index]  unmapped group id (<query_label>.iteration_XXXXXXXX.unmapped_j) => member ids
        """
        ids = query_ids
        sizes = ones(Int,length(ids)) # leaves

        query_mappings_list = Dict{String,String}[]
        unmapped_groups_list = Dict{String,Vector{String}}[]
        previous_level = nothing # the level below; only needed to compute levels without binaries
        if !isfile(reference_pathlist[1])
            leaves = [load_milk_binary_matrix(path) for path in reference_leaves_pathlist]
            leaf_ids = reduce(vcat,first.(leaves))
            previous_level = ReferenceLevel("leaves",leaf_ids,reduce(hcat,last.(leaves)),ones(Int,length(leaf_ids)),zeros(Int,length(leaf_ids)),0.0)
        end
        for (index,reference_path) in enumerate(reference_pathlist)
            iteration = index-1 # MILK iterations are 0-based
            reference = isfile(reference_path) ? load_reference_level(reference_path) : compute_reference_level(reference_path,previous_level)
            previous_level = reference
            mappings,unmapped_groups = place_query_objects_to_reference_level(reference,X_query,sizes,rng;partition_size=partition_size)

            query_mappings = Dict{String,String}()
            for (i,id) in enumerate(ids)
                if mappings[i] > 0
                    query_mappings[id] = reference.ids[mappings[i]]
                end
            end

            unmapped_ids = ["$(query_label).iteration_$(lpad(iteration,8,'0')).unmapped_$(j)" for j in eachindex(unmapped_groups)]
            unmapped_members = Dict(unmapped_id => ids[group] for (unmapped_id,group) in zip(unmapped_ids,unmapped_groups))

            if iteration == length(reference_pathlist)-1 # top iteration: whatever is left joins the root
                for members in values(unmapped_members), member in members
                    query_mappings[member] = reference.ids[1]
                end
                unmapped_groups = Vector{Int}[]
                unmapped_members = Dict{String,Vector{String}}()
            end
            push!(query_mappings_list,query_mappings)
            push!(unmapped_groups_list,unmapped_members)
            @info "\tIteration $(iteration): $(length(ids)) query objects; $(length(query_mappings)) mapped to the reference; $(length(unmapped_groups)) unmapped groups"

            if isempty(unmapped_groups)
                break
            end
            unmapped_sizes = [sum(sizes[unmapped_group]) for unmapped_group in unmapped_groups]
            X_query = stack(_weighted_centroid(X_query[:,unmapped_group],sizes[unmapped_group]) for unmapped_group in unmapped_groups)
            sizes = unmapped_sizes
            ids = unmapped_ids
        end
        while length(query_mappings_list) < length(reference_pathlist) # placement stopped early: nothing changes above
            push!(query_mappings_list,Dict{String,String}())
            push!(unmapped_groups_list,Dict{String,Vector{String}}())
        end
        return query_mappings_list,unmapped_groups_list
    end

    function recompute_backbone(reference_pathlist,query_mappings_list,unmapped_groups_list,output_dir,label,stage,ids_pathlist)
        """
        New tree levels, bottom-up: each reference node's members plus the query objects mapped to it, and the unmapped
        groups as new nodes. write_tree_level recomputes centroid, size, spread and radius (existing nodes keep their stage,
        new nodes get `stage`). ids_pathlist: input binaries of every cell placed so far.
        """
        tmp_ids_pathlist = String[]
        previous_ids_pathlist = ids_pathlist
        for (index,reference_path) in enumerate(reference_pathlist)
            iteration = index-1 # MILK iterations are 0-based
            groups = load_groups_as_dictionary(jsonl_path(reference_path))["groups"]
            for (query_id,reference_id) in query_mappings_list[index]
                push!(groups[reference_id],query_id)
            end
            merge!(groups,unmapped_groups_list[index])
            node_stage_dict = isfile(reference_path) ? load_milk_group_information(reference_path).stage : Dict{String,Int}()

            # written as .tmp, moved once every iteration is written (old levels stay intact if this fails midway)
            tmp_ids_path = replace(tree_level_path(output_dir,label,iteration),r"\.ids$" => ".tmp.ids")
            write_tree_level(previous_ids_pathlist,groups,get_threshold(reference_path),tmp_ids_path,stage;node_stage_dict=node_stage_dict)
            push!(tmp_ids_pathlist,tmp_ids_path)
            previous_ids_pathlist = [tmp_ids_path]
        end

        for (index,tmp_ids_path) in enumerate(tmp_ids_pathlist)
            ids_path = tree_level_path(output_dir,label,index-1)
            mv_milk_files(tmp_ids_path,ids_path;force=true)
            mv(jsonl_path(tmp_ids_path),jsonl_path(ids_path);force=true)
            mv(metadata_path(jsonl_path(tmp_ids_path)),metadata_path(jsonl_path(ids_path));force=true)
        end
        @info "Recomputed: $(length(reference_pathlist)) iterations written to $(output_dir)"
    end

    function write_query_metadata(path,query_mappings_list,unmapped_groups_list,query_ids)
        placement_dict = Dict{String,Tuple{Int,String,String}}()
        for index in length(query_mappings_list):-1:1 # top-down
            for (query_id,reference_id) in query_mappings_list[index]
                placement_dict[query_id] = (index-1,reference_id,"")
            end
            for (unmapped_id,members) in unmapped_groups_list[index]
                iteration,reference_id,top_id = placement_dict[unmapped_id]

                # top unmapped group of its query-only branch (empty if the object was mapped directly)
                branch_id = isempty(top_id) ? unmapped_id : top_id
                for member in members
                    placement_dict[member] = (iteration,reference_id,branch_id)
                end
            end
        end
        open(path,"w") do io
            for query_id in query_ids
                # julia splat -> [query_id,iteration,reference_id, branch_id]
                println(io,join([query_id,placement_dict[query_id]...],'\t'))
            end
        end
    end

    function placement(reference_dir,output_dir,label,query_ids_path,reference_leaves_pathlist;query_label="query",stage=0,seed=0,partition_size=10000)
        reference_pathlist = get_reference_pathlist(reference_dir,label)
        query_ids,X_query = load_milk_binary_matrix(query_ids_path)
        mkpath(output_dir)
        rng = Xoshiro(seed)

        query_mappings_list,unmapped_groups_list = place_query_objects_to_reference(reference_pathlist,reference_leaves_pathlist,query_ids,X_query,query_label,rng;partition_size=partition_size)
        placement_metadata_path = joinpath(output_dir,"$(label).$(query_label).placement.tsv")
        write_query_metadata(placement_metadata_path,query_mappings_list,unmapped_groups_list,query_ids)
        recompute_backbone(reference_pathlist,query_mappings_list,unmapped_groups_list,output_dir,label,stage,vcat(reference_leaves_pathlist,[query_ids_path]))
        return query_mappings_list,unmapped_groups_list
    end
end