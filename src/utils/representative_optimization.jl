module RepresentativeOptimization

    using Statistics

    export optimize_representatives

    function compute_centroid(medoid_dict)
        arr = reduce(hcat,values(medoid_dict))
        return vec(mean(arr,dims=2))
    end

    function get_medoid_dictionary(group,representative_dict,cache_dict)
        medoid_dict = Dict{String,Vector{Float32}}() # this is for a single group
        cache_set = Set{String}()
        for sample_id in group
            if haskey(representative_dict,sample_id)
                medoid_dict[sample_id] = representative_dict[sample_id]
            elseif !isnothing(cache_dict) && haskey(cache_dict,sample_id)
                medoid_dict[sample_id] = cache_dict[sample_id]
                push!(cache_set,sample_id)
            end
        end
        return medoid_dict,cache_set
    end

    function determine_medoid(medoid_dict,centroid_vec,cache_set,distance_function)
        """
        The optimized medoid cannot be the cached set of object, it must be an
        object that's present in the current recursive iteration.
        """
        medoid_id = nothing
        medoid_vec = nothing
        medoid_dist = Inf
        n_comparisons = 0
        for (sample_id,vec) in medoid_dict
            if sample_id in cache_set continue end
            dist = distance_function(centroid_vec,vec)
            n_comparisons += 1
            if dist < medoid_dist
                medoid_id   = sample_id
                medoid_vec  = vec
                medoid_dist = dist
            end
        end
        return medoid_id,medoid_vec,n_comparisons
    end

    function optimize_representatives(data_dict,groups,cache_dict,distance_function)
        """
        TODO: determine if there are enough cases where medoid_dict (and associated cached objects) are 
        <= 2, as that's the trivial case. May get marginal optimization by avoiding those computing tasks.
        """
        optimized_dict = Dict{String,Vector{Float32}}()
        optimization_set = Set{String}()
        n_comparisons = 0
        for (representative_id,group) in groups
            if length(group) == 1
                optimized_dict[representative_id] = data_dict[representative_id]
            else
                medoid_dict,cache_set = get_medoid_dictionary(group,data_dict,cache_dict)
                centroid_vec = compute_centroid(medoid_dict)
                medoid_id,medoid_vec,n_centroid_comps = determine_medoid(medoid_dict,centroid_vec,cache_set,distance_function)
                n_comparisons += n_centroid_comps
                optimized_dict[medoid_id] = medoid_vec
                push!(optimization_set,medoid_id)
            end
        end
        return optimized_dict,optimization_set,n_comparisons
    end

end