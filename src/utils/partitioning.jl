module Partitioning

    using Glob
    using Random
    using Logging
    using StatsBase
    using Distributed
    using Statistics
    using Distances
    using LinearAlgebra

    using ..FileHandling: binary_path,cache_meta_path,rm_milk_files,count_rows,
                          open_file_limit,read_binary_rows,open_milk_files,close_milk_files

    export partition_input_file

    # k-means distances: squared Euclidean gives the same nearest centroid as Euclidean and is cheaper
    const METRIC = Dict(
        "euclidean"   => SqEuclidean(),
        "hamming"     => SqEuclidean(),
        "jaccard"     => SqEuclidean(),
        "manhattan"   => Cityblock(),
        "cosine"      => CosineDist(),
        "correlation" => CorrDist()
    )

    mutable struct Node
        centroid::Vector{Float32}
        X_children::Matrix{Float32} # d x b child centroids (empty for leaves)
        children::Vector{Int}
    end

    function _is_leaf(node)
        return isempty(node.children)
    end

    function _partition_path(partition_dir,label,p)
        return joinpath(partition_dir,"$(label).partition_$(lpad(p,8,'0')).ids")
    end

    function sample_objects(ids_path,n_input_objects,n_input_rows,partition_size,rng;
                            sample_bytes_cutoff=1<<30,per_centroid_coverage,
                            min_sample_size=200000,min_chunk_size=25000)

        bin_path = binary_path(ids_path)
        m = div(filesize(bin_path),4*n_input_objects) # dimensionality

        sample_size = min(
            n_input_objects,
            div(sample_bytes_cutoff,4*m), # memory-bound to 1GB of Float32 values
            max(
                per_centroid_coverage*cld(n_input_rows,partition_size), # 100 objects per centroid (partition)
                min_sample_size 
            )
        )

        sampled_indices = sort!(sample(rng,1:n_input_objects,sample_size;replace=false))
        n_chunks = min(nworkers(),cld(sample_size,min_chunk_size))
        chunks   = collect(Iterators.partition(sampled_indices,cld(sample_size,n_chunks)))
        X_chunks = pmap(rows -> read_binary_rows(bin_path,rows,m),chunks)
        return reduce(hcat,X_chunks) # m x sample_size, column-major
    end

    # ---------------------------------------------------------------- k-means

    # sample rows for cosine/correlation: center (correlation) and unit-normalize
    function normalize_embeddings!(X,metric)
        if !(metric in ("cosine","correlation"))
            return X
        end
        @inbounds @views for j in axes(X,2)
            vec = X[:,j]
            if metric == "correlation"
                vec .-= mean(vec)
            end
            norm_val = norm(vec)
            if norm_val > 0
                vec ./= norm_val
            end
        end
        return X
    end

    # NaN/Inf (e.g. a zero vector under cosine) -> Inf32; tiny negatives from rounding -> 0
    function _clean_cost(x)
        return isfinite(x) ? max(x,0) : Inf32
    end

    # k-means++ if distance_function == squared Euclidean
    function _centroid_initialization(X,k,distance_function,rng)
        m,n = size(X)
        C = Matrix{Float32}(undef,m,k)

        C[:,1] = X[:,rand(rng,1:n)]
        costs = _clean_cost.(colwise(distance_function,X,C[:,1]))

        for j in 2:k
            sampling_costs = ifelse.(isfinite.(costs),costs,0)
            cumulative_costs = cumsum(sampling_costs)
            total_cost = cumulative_costs[end]
            if isfinite(total_cost) && total_cost > 0
                i = searchsortedfirst(cumulative_costs,(1-rand(rng))*total_cost)
            else
                # every object sits on a centroid (or invalid costs): sample uniformly
                i = rand(rng,1:n)
            end
            C[:,j] = X[:,i]
            costs .= min.(costs,_clean_cost.(colwise(distance_function,X,C[:,j])))
        end
        return C
    end

    function _assignment_update!(assignments,C,X,distance_function)
        D = pairwise(distance_function,C,X;dims=2) # k x n
        converged = true
        @inbounds for i in axes(X,2)
            j = argmin(view(D,:,i))
            if j != assignments[i]
                converged = false
                assignments[i] = j
            end
        end
        return converged
    end

    function _centroid_update!(C,X,assignments,rng)
        m,n = size(X)
        k = size(C,2)
        X_sums = zeros(Float32,m,k)
        sizes = zeros(Int,k)
        @inbounds for i in 1:n
            j = assignments[i]
            @views X_sums[:,j] .+= X[:,i]
            sizes[j] += 1
        end
        for j in 1:k
            if sizes[j] > 0
                @views C[:,j] .= X_sums[:,j] ./ sizes[j] # arithmetic mean
            else
                C[:,j] .= X[:,rand(rng,1:n)]
            end
        end
    end

    function generalized_kmeans(X,k,distance_function,rng;max_iter=20)
        C = _centroid_initialization(X,k,distance_function,rng)
        assignments = zeros(Int,size(X,2))
        for iter in 1:max_iter
            converged = _assignment_update!(assignments,C,X,distance_function)
            if iter == max_iter || converged
                break
            end
            _centroid_update!(C,X,assignments,rng)
        end
        return C,assignments
    end

    function _get_branching_factor(estimated_size,partition_size,n_samples;branching_factor=10)
        return min(
            branching_factor,
            ceil(Int,estimated_size/partition_size),
            n_samples
        )
    end

    function build_initial_tree(X,n_input_rows,partition_size,max_leaves,distance_function,rng;oversized_tolerance)
        sample_size = size(X,2)

        root = Node(Float32[],zeros(Float32,size(X,1),0),Int[]) # centroid, X_children, children
        nodes = [root]
        queue = [(1,collect(1:sample_size))]
        n_leaves = 1
        while !isempty(queue)
            node_id,indices = popfirst!(queue)
            node = nodes[node_id]
            estimated_size = length(indices)*n_input_rows/sample_size
            k = _get_branching_factor(estimated_size,partition_size,length(indices))
            if estimated_size <= oversized_tolerance*partition_size || n_leaves - 1 + k > max_leaves || k < 2
                continue
            end

            C,assignments = generalized_kmeans(X[:,indices],k,distance_function,rng)
            viable_assignments = sort(unique(assignments))
            if length(viable_assignments) < 2
                continue # no separation (e.g. duplicates)
            end

            node.X_children = C[:,viable_assignments]
            for j in viable_assignments
                child_indices = indices[assignments .== j]
                child_node = Node(C[:,j],zeros(Float32,size(X,1),0),Int[]) # centroid, X_children, children
                push!(nodes,child_node)
                push!(node.children,length(nodes))
                push!(queue,(length(nodes),child_indices))
            end
            n_leaves += length(viable_assignments) - 1
        end
        return nodes
    end

    function get_leaf_neighbors(nodes,distance_function;n_neighbors=10)
        ids = findall(_is_leaf,nodes)
        n_leaves = length(ids)
        m = size(nodes[1].X_children,1) # dimensionality
        X_leaves = zeros(Float32,m,n_leaves)
        leaf_positions = zeros(Int,length(nodes)) # node id -> leaf position (0 for internal nodes)
        for (j,id) in enumerate(ids)
            if !isempty(nodes[id].centroid)
                X_leaves[:,j] = nodes[id].centroid
            end
            leaf_positions[id] = j
        end

        k = min(n_neighbors,n_leaves-1)
        neighbors_dict = Dict{Int,Vector{Int}}()
        for block in Iterators.partition(1:n_leaves,1024)
            D = pairwise(distance_function,X_leaves,X_leaves[:,block],dims=2)
            for (b,j) in enumerate(block)
                D[j,b] = -Inf # in event of ties, force self to be first
                neighbors_dict[j] = partialsortperm(view(D,:,b),1:k+1)
            end
        end
        return X_leaves,leaf_positions,neighbors_dict
    end

    # ids grouped by their current assignment
    function _group_by(ids,assignments)
        assignments_dict = Dict{Int,Vector{Int}}()
        for id in ids
            assignment = assignments[id]
            if !haskey(assignments_dict,assignment)
                assignments_dict[assignment] = Int[]
            end
            push!(assignments_dict[assignment],id)
        end
        return assignments_dict
    end


    function greedy_hierarchical_mapping(X,nodes,distance_function)
        n = size(X,2)
        assignments = zeros(Int,n) # leaf node id per object
        objects_at = Dict(1 => collect(1:n)) # node id -> objects currently there
        for node_id in eachindex(nodes)
            objects = pop!(objects_at,node_id,Int[])
            if isempty(objects)
                continue
            end
            node = nodes[node_id]
            if _is_leaf(node)
                assignments[objects] .= node_id
                continue
            end
            D = pairwise(distance_function,node.X_children,X[:,objects],dims=2)
            for (k,i) in enumerate(objects)
                child_id = node.children[argmin(view(D,:,k))]
                push!(get!(objects_at,child_id,Int[]),i)
            end
        end
        return assignments
    end

    # greedy leaf, then the nearest among that leaf's neighbouring leaves
    function mapping(X,nodes,X_leaves,leaf_positions,neighbors_dict,distance_function)
        node_ids = greedy_hierarchical_mapping(X,nodes,distance_function)
        assignments = [leaf_positions[id] for id in node_ids]
        for (j,objects) in _group_by(1:size(X,2),assignments)
            candidates = neighbors_dict[j]
            D = pairwise(distance_function,X_leaves[:,candidates],X[:,objects],dims=2)
            for (k,i) in enumerate(objects)
                assignments[i] = candidates[argmin(view(D,:,k))]
            end
        end
        return assignments
    end

    # objects (columns of the batch) appended to one set of MILK files, with their cached cells
    function _write_objects(files,objects,ids,X,groups,cached_metadata,cached_binaries,cache_flag)
        print(files.ids,join(ids[objects],'\n'),'\n')
        write(files.bin,X[:,objects])
        print(files.group,join(groups[objects],'\n'),'\n')
        if cache_flag
            for j in objects
                for line in cached_metadata[j]
                    println(files.cache_meta,line)
                end
                write(files.cache_bin,cached_binaries[j])
            end
        end
    end

    function write_partitions(ids_path,n_input_objects,m,partition_pathlist,threshold_pathlist,threshold_index_map,
                              nodes,X_leaves,leaf_positions,neighbors_dict,distance_function;chunk_cutoff_bytes=1<<26)
        n_partitions = length(partition_pathlist)
        partition_sizes = zeros(Int,n_partitions)
        cache_flag = isfile(cache_meta_path(ids_path))

        inputs = open_milk_files(ids_path,"r";cache=cache_flag)
        partitioned_inputs = [open_milk_files(path,"w";cache=cache_flag) for path in partition_pathlist]
        threshold_inputs = [open_milk_files(path,"w";cache=cache_flag) for path in threshold_pathlist]
        
        chunk_size = max(256,div(chunk_cutoff_bytes,4*m)) # one chunk per worker per batch
        batch_size = chunk_size*nworkers()
        pool = CachingPool(workers()) # sends the tree to each worker once
        route_chunk = X_chunk -> mapping(X_chunk,nodes,X_leaves,leaf_positions,neighbors_dict,distance_function)

        for batch in Iterators.partition(1:n_input_objects,batch_size)
            n_batch = length(batch)
            X = read!(inputs.bin,Matrix{Float32}(undef,m,n_batch))
            ids = [readline(inputs.ids) for _ in 1:n_batch]
            groups = [readline(inputs.group) for _ in 1:n_batch]
            for j in 1:n_batch
                if !startswith(groups[j],ids[j]*"\t")
                    error("Row mismatch between .ids and .group.tsv: $(ids[j]) vs $(groups[j])")
                end
            end

            n_object_rows = ones(Int,n_batch)
            cached_metadata = Vector{Vector{String}}(undef,n_batch)
            cached_binaries = Vector{Vector{Float32}}(undef,n_batch)
            if cache_flag
                for j in 1:n_batch
                    n_object_rows[j] = parse(Int,split(groups[j],'\t')[4])
                    cached_metadata[j] = [readline(inputs.cache_meta) for _ in 1:n_object_rows[j]]
                    cached_binaries[j] = read!(inputs.cache_bin,Vector{Float32}(undef,m*n_object_rows[j]))
                end
            end

            # batch routing: one chunk of columns per worker. Raw vectors are routed: cosine and correlation
            # distances ignore scale, so only the sample used to build the tree needed normalizing
            chunks = collect(Iterators.partition(1:n_batch,chunk_size))
            assignments = reduce(vcat,pmap(route_chunk,pool,[X[:,chunk] for chunk in chunks]))

            for (p,objects) in _group_by(1:n_batch,assignments)
                _write_objects(partitioned_inputs[p],objects,ids,X,groups,cached_metadata,cached_binaries,cache_flag)
                partition_sizes[p] += sum(n_object_rows[objects])
            end

            threshold_assignments = [get(threshold_index_map,index,0) for index in batch]
            for (t,objects) in _group_by(findall(>(0),threshold_assignments),threshold_assignments)
                _write_objects(threshold_inputs[t],objects,ids,X,groups,cached_metadata,cached_binaries,cache_flag)
            end
        end

        clear!(pool)
        close_milk_files(inputs)
        foreach(close_milk_files,partitioned_inputs)
        foreach(close_milk_files,threshold_inputs)
        return partition_sizes
    end

    function partition_file!(partition_counter,ids_path,partition_dir,label,invariant_args,rng;
                             per_centroid_coverage=100,n_reserved_files=64,oversized_tolerance=1.25,
                             min_shrinkage_factor=0.95,threshold_dir=nothing)
        start_time = time()
        partition_size = invariant_args["partition-size"]
        metric = invariant_args["metric"]
        distance_function = METRIC[metric]

        n_input_objects = countlines(ids_path)
        n_input_rows = count_rows(ids_path) # includes cached objects (equal to n_input_objects if no cache)
        files_per_partition = isfile(cache_meta_path(ids_path)) ? 5 : 3

        threshold_pathlist = String[]
        threshold_index_map = Dict{Int,Int}() # input object index -> threshold file
        if !isnothing(threshold_dir) # top-level run only
            n_threshold_files = min(invariant_args["batch-size"],cld(n_input_rows,partition_size))
            n_threshold_objects = min(n_input_objects,round(Int,n_threshold_files*partition_size*n_input_objects/n_input_rows))
            for (k,index) in enumerate(sample(rng,1:n_input_objects,n_threshold_objects;replace=false))
                threshold_index_map[index] = mod1(k,n_threshold_files)
            end
            for t in 1:n_threshold_files
                push!(threshold_pathlist,_partition_path(threshold_dir,label,t))
            end
        end

        X_sample = sample_objects(ids_path,n_input_objects,n_input_rows,partition_size,rng;per_centroid_coverage=per_centroid_coverage)
        m = size(X_sample,1)
        normalize_embeddings!(X_sample,metric)

        max_leaves = max(
            2,
            min(
                div(size(X_sample,2),per_centroid_coverage),
                div(open_file_limit()-(files_per_partition*length(threshold_pathlist)+n_reserved_files),files_per_partition)
            )
        )
        nodes = build_initial_tree(X_sample,n_input_rows,partition_size,max_leaves,distance_function,
                                   rng;oversized_tolerance=oversized_tolerance)
        X_leaves,leaf_positions,neighbors_dict = get_leaf_neighbors(nodes,distance_function)
        tree_time = round(time()-start_time,digits=1)

        n_partitions = size(X_leaves,2)
        partition_pathlist = String[]
        for _ in 1:n_partitions
            partition_counter[] += 1
            push!(partition_pathlist,_partition_path(partition_dir,label,partition_counter[]))
        end
        partition_sizes = write_partitions(ids_path,n_input_objects,m,partition_pathlist,threshold_pathlist,threshold_index_map,
                                           nodes,X_leaves,leaf_positions,neighbors_dict,distance_function)
        write_time = round(time()-start_time-tree_time,digits=1)

        @info "\t\t$(basename(ids_path)): $(n_input_rows) rows -> $(count(>(0),partition_sizes)) partitions (largest $(maximum(partition_sizes))); tree $(tree_time)s, write $(write_time)s"
        flush(stdout)

        for (partition_path,partition_rows) in zip(partition_pathlist,partition_sizes)
            if partition_rows == 0
                rm_milk_files(partition_path) # empty leaf
            elseif partition_rows <= oversized_tolerance*partition_size
                continue
            elseif partition_rows <= min_shrinkage_factor*n_input_rows
                partition_file!(partition_counter,partition_path,partition_dir,label,invariant_args,rng)
                rm_milk_files(partition_path)
            else
                @warn "Could not split $(basename(ids_path)): a partition kept $(partition_rows) of $(n_input_rows) rows; keeping it as one oversized partition"
            end
        end
    end

    function partition_input_file(input_path,label,invariant_args)
        partition_dir = joinpath(invariant_args["output-dir"],"$(label).split")
        threshold_dir = joinpath(partition_dir,"threshold_sample")
        mkpath(threshold_dir)

        rng = Xoshiro(invariant_args["seed"])
        partition_counter = Ref(0)
        partition_file!(partition_counter,input_path,partition_dir,label,invariant_args,rng;threshold_dir=threshold_dir)
        threshold_pathlist = sort(glob("*.ids",threshold_dir))
        return partition_dir,threshold_pathlist
    end

end