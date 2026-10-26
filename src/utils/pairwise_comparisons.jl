module PairwiseComparisons

    using StatsBase
    using Distances
    using Distributed

    using ..FileHandling: load_milk_binaries
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

    function compute_pairwise_distance_matrix(data_dict,p,distance_function,compute_threshold=true)
        """
        When working with lower-dimensional representations of gene expression space,
        distribution of distance matrix computation across files with single CPU distance matrix 
        computation is much faster than sequentially computing distance matrices with multiple
        threads.
        'p' denotes percentile.
        """
        ids = sort!(collect(keys(data_dict)))
        index_map = Dict( ids[i] => i for i in eachindex(ids) )
        n = length(ids)
        X = stack(data_dict[id] for id in ids) # d × n, one column per object
        metric = get(PAIRWISE_METRIC,distance_function,nothing)
        D = isnothing(metric) ?
            Float32[distance_function(view(X,:,i),view(X,:,j)) for i in 1:n, j in 1:n] :
            pairwise(metric,X;dims=2)
        D = eltype(D) === Float32 ? D : Float32.(D)
        for i in 1:n
            D[i,i] = 0
        end
        threshold = compute_threshold ? percentile([D[i,j] for j in 2:n for i in 1:j-1],p) : nothing
        n_comparisons = div(n*(n-1),2)
        return D,index_map,threshold,n_comparisons
    end

    function determine_percentile_threshold_process(input_path,metric,p,linkage)
        data_dict,groupsize_dict = load_milk_binaries(input_path)
        if length(data_dict) < 2
            return nothing
        end
        distance_function = map_distance_function(metric)
        # D,index_map,threshold,_ = compute_pairwise_distance_matrix(data_dict,p,distance_function)
        D,index_map,threshold,_ = compute_pairwise_distance_matrix(data_dict,p,distance_function,!startswith(linkage,"ward"))
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