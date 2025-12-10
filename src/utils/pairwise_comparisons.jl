module PairwiseComparisons

    using StatsBase
    using Distances
    using Distributed

    using ..FileHandling: load_input_array_as_dictionary

    export map_distance_function,
           compute_pairwise_distance_matrix,
           determine_percentile_threshold

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

    function compute_pairwise_distance_matrix(data_dict,p,distance_function)
        """
        When working with lower-dimensional representations of gene expression space,
        distribution of distance matrix computation across files with single CPU distance matrix 
        computation is much faster than sequentially computing distance matrices with multiple
        threads.
        'p' denotes percentile.
        """
        ids = collect(keys(data_dict))
        index_map = Dict( ids[i] => i for i in eachindex(ids) )
        n = length(ids)
        D = zeros(Float32,n,n)
        indices = [ (ids[i],ids[j]) for i in 1:n-1 for j in i+1:n ]
        L = length(indices)
        distances = Vector{Float32}(undef,L)
        @inbounds @simd for x in eachindex(indices)
            id_i,id_j = indices[x]
            i = index_map[id_i]
            j = index_map[id_j]
            d = distance_function(data_dict[id_i],data_dict[id_j])
            D[i,j] = d
            D[j,i] = d
            distances[x] = d
        end
        threshold = percentile(distances,p)
        n_comparisons = length(distances)
        return D,index_map,threshold,n_comparisons
    end

    function determine_percentile_threshold_process(input_path,metric,p)
        data_dict = load_input_array_as_dictionary(input_path)
        distance_function = map_distance_function(metric)
        _,_,threshold,_ = compute_pairwise_distance_matrix(data_dict,p,distance_function)
        return threshold
    end

    function determine_percentile_threshold(pathlist,invariant_args)
        """ distributed
        """
        metric = invariant_args["metric"]
        p = invariant_args["percentile"]
        thresholds = pmap(path -> determine_percentile_threshold_process(path,metric,p),pathlist)
        GC.gc() 
        return minimum(thresholds)
    end

    @inline function update_closest_representative(dist,threshold,representative_id,closest_id,closest_dist,specificity)
        if dist <= threshold
            specificity += 1
            if dist < closest_dist
                return representative_id,dist,specificity
            end
        end
        return closest_id,closest_dist,specificity
    end

    function map_object_to_representative(vec,representative_dict,threshold,distance_function)
        """
        After representative optimization, re-map objects.
        Map each candidate (denoted by its index, i) to the representative it is closest to
        Returns nothing if there are no representatives it is similar to, as defined by
        the threshold. Computing pairwise distances.
        """
        closest_id = nothing
        closest_dist = Inf
        specificity = 0
        n_comparisons = 0
        for (representative_id,representative_vec) in representative_dict
            dist = distance_function(vec,representative_vec)
            n_comparisons += 1
            closest_id,closest_dist,specificity = update_closest_representative(dist,threshold,representative_id,
                                                                                closest_id,closest_dist,specificity)
        end
        return closest_id,closest_dist,n_comparisons,specificity
    end

    function map_object_to_representative_precomputed(i,representative_dict,D,index_map,threshold)
        """
        After representative optimization, re-map objects.
        Map each candidate (denoted by its index, i) to the representative it is closest to
        Returns nothing if there are no representatives it is similar to, as defined by
        the threshold. Using precomputed all-all distance matrix.
        """
        closest_id = nothing
        closest_dist = Inf
        specificity = 0
        for representative_id in keys(representative_dict)
            j = index_map[representative_id]
            @inbounds dist = D[i,j]
            closest_id,closest_dist,specificity = update_closest_representative(dist,threshold,representative_id,
                                                                                closest_id,closest_dist,specificity)
        end
        return closest_id,closest_dist,specificity
    end

end