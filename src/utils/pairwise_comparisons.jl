module PairwiseComparisons

    using StatsBase
    using Distances
    using Distributed

    using ..FileHandling: load_milk_binaries,load_milk_cache,cache_meta_path
    using ..Ward: compute_ward_cost_matrix,ward_percentile_threshold

    export map_distance_function,
           compute_pairwise_distance_matrix,
           determine_percentile_threshold,
           compute_euclidean_distance

    function compute_cosine_distance(vec_i,vec_j)
        return cosine_dist(vec_i,vec_j)
    end

    function compute_euclidean_distance(vec_i,vec_j)
        return euclidean(vec_i,vec_j)
    end

    function compute_manhattan_distance(vec_i,vec_j)
        return cityblock(vec_i,vec_j)
    end

    function compute_hamming_distance(vec_i,vec_j)
        return hamming(vec_i,vec_j)
    end

    function compute_jaccard_distance(vec_i,vec_j)
        return jaccard(vec_i,vec_j)
    end

    function compute_correlation_distance(vec_i,vec_j)
        return corr_dist(vec_i,vec_j)
    end

    function map_distance_function(metric)
        if metric == "cosine"
            distance_function = compute_cosine_distance
        elseif metric == "euclidean"
            distance_function = compute_euclidean_distance
        elseif metric == "manhattan"
            distance_function = compute_manhattan_distance
        elseif metric == "hamming"
            distance_function = compute_hamming_distance
        elseif metric == "jaccard"
            distance_function = compute_jaccard_distance
        elseif metric == "correlation"
            distance_function = compute_correlation_distance
        else
            throw("Specified metric ($metric) not supported!")
        end
        return distance_function
    end

    const PAIRWISE_METRIC = Dict(
        compute_euclidean_distance => Euclidean(),
        compute_cosine_distance => CosineDist(),
        compute_manhattan_distance => Cityblock(),
        compute_hamming_distance => Hamming(),
        compute_jaccard_distance => Jaccard(),
        compute_correlation_distance => CorrDist(),
    )

    function _distance_matrix(X,distance_function)
        # columns of X are vectors
        metric = get(PAIRWISE_METRIC,distance_function,nothing)
        n = size(X,2)
        D = isnothing(metric) ?
            Float32[distance_function(view(X,:,i),view(X,:,j)) for i in 1:n, j in 1:n] :
            pairwise(metric,X;dims=2)
        return eltype(D) === Float32 ? D : Float32.(D)
    end

    function _cached_distance_matrix(ids,cache_dict,distance_function)
        """
        Group-by-group distance matrix computed from each group's representative and its cached cells.

        1. Compute all-vs-all distances over every group's representative and cached cells.
        2. Collapse to group-by-group: entry (g,h) is the weighted average of the distances between
        g's representative + cached cells and h's representative + cached cells. Each pair is weighted
        by the fraction of g the first cell stands for times the fraction of h the second cell stands for.
        A cell's fraction is its subgroup_size divided by its group's total size, so each group's
        fractions sum to 1.

        With cache_size = 0 each group contributes only its representative (fraction 1), so entry (g,h)
        is the plain representative-to-representative distance.
        """
        vecs = Vector{Float32}[]
        groups = Int[] # index (in ids) of the group each cell belongs to
        weights = Float64[] # fraction of its group the cell stands for
        for (g,id) in enumerate(ids)
            cache = cache_dict[id]
            group_size = sum(cache.subgroup_sizes)
            for k in eachindex(cache.cell_ids)
                push!(vecs,cache.vecs[k])
                push!(groups,g)
                push!(weights,cache.subgroup_sizes[k]/group_size)
            end
        end
        n_cells = length(vecs)
        n_groups = length(ids)
    
        D = _distance_matrix(stack(vecs),distance_function)

        # each cell's weighted average distance to each group
        cell2group_dist = zeros(Float64,n_cells,n_groups)
        for b in 1:n_cells
            cell2group_dist[:,groups[b]] .+= weights[b] .* view(D,:,b)
        end

        # each group's weighted average distance to each group
        group_dist = zeros(Float64,n_groups,n_groups)
        for a in 1:n_cells
            group_dist[groups[a],:] .+= weights[a] .* view(cell2group_dist,a,:)
        end
    
        # same average summed in two orders can differ in the last digit; average both halves so D[g,h] == D[h,g]
        group_dist = (group_dist .+ transpose(group_dist)) ./ 2
        return Float32.(group_dist),n_cells
    end

    function compute_pairwise_distance_matrix(data_dict,p,distance_function,compute_threshold=true;cache_dict=nothing)
        """
        'p' denotes percentile. With cache_dict, distances are averaged over each group's cached cells;
        without it, distances are between representatives.
        """
        ids = sort!(collect(keys(data_dict)))
        index_map = Dict( ids[i] => i for i in eachindex(ids) )
        n = length(ids)
        if isnothing(cache_dict)
            D = _distance_matrix(stack(data_dict[id] for id in ids),distance_function)
            n_rows = n
        else
            D,n_rows = _cached_distance_matrix(ids,cache_dict,distance_function)
        end
        for i in 1:n
            D[i,i] = 0
        end
        threshold = compute_threshold ? percentile([D[i,j] for j in 2:n for i in 1:j-1],p) : nothing
        n_comparisons = div(n_rows*(n_rows-1),2)
        return D,index_map,threshold,n_comparisons
    end

    function determine_percentile_threshold_process(input_path,metric,p,linkage)
        data_dict,groupsize_dict = load_milk_binaries(input_path)
        if length(data_dict) < 2
            return nothing
        end
        distance_function = map_distance_function(metric)
        if (!startswith(linkage,"ward") && isfile(cache_meta_path(input_path)))
            cache_dict = load_milk_cache(input_path,data_dict,groupsize_dict)
        else
            cache_dict = nothing
        end
        D,index_map,threshold,_ = compute_pairwise_distance_matrix(data_dict,p,distance_function,!startswith(linkage,"ward");cache_dict=cache_dict)
        if startswith(linkage,"ward")
            balanced = linkage == "ward_balanced"
            ids = Vector{String}(undef,length(index_map))
            for (id,i) in index_map
                ids[i] = id
            end
            sizes = [groupsize_dict[id] for id in ids]
            cost_sizes = balanced ? ones(Int,length(sizes)) : copy(sizes)
            threshold = ward_percentile_threshold(compute_ward_cost_matrix(D,cost_sizes),p)
        end
        return threshold
    end

    function determine_percentile_threshold(pathlist,invariant_args)
        """ distributed
        """
        metric = invariant_args["metric"]
        p = invariant_args["percentile"]
        linkage = invariant_args["linkage"]
        thresholds = pmap(path -> determine_percentile_threshold_process(path,metric,p,linkage),pathlist)
        filtered_thresholds = filter(!isnothing,thresholds)
        if isempty(filtered_thresholds)
            error("Unable to compute threshold due to all partitions having too few objects.")
        end
        GC.gc() 
        return minimum(filtered_thresholds)
    end

end