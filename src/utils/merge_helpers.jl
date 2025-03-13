
include("pairwise_comparisons.jl")

function stratification(data_dict,threshold,distance_function)
    representatives_dict = Dict{String,Vector{Float32}}()
    groups = Dict{String,Vector{String}}()
    x = 0
    for (candidate_id,candidate_vec) in data_dict
        closest_id = nothing
        closest_dist = Inf
        for (representative_id,representative_vec) in representatives_dict
            dist = distance_function(candidate_vec,representative_vec)
            x += 1
            if (dist <= threshold) && (dist < closest_dist)
                closest_id = representative_id
                closest_dist = dist
            end
        end
        if isnothing(closest_id)
            representatives_dict[candidate_id] = candidate_vec
            groups[candidate_id] = [candidate_id]
        else
            push!(groups[closest_id],candidate_id)
        end
    end
    return representatives_dict,groups,x
end

function stratification_predefined_medoids(;data_dict,representative_dict,threshold,distance_function)
    groups = Dict{String,Vector{String}}()
    unmapped_dict = Dict{String,Vector{Float32}}()
    distances_dict = Dict{String,Float32}()
    specificity_dict = Dict{String,Int32}()
    x = 0
    for (candidate_id,candidate_vec) in data_dict
        if haskey(representative_dict,candidate_id) continue end
        closest_id = nothing
        closest_dist = Inf
        specificity = 0
        for (representative_id,representative_vec) in representative_dict
            get!(groups,representative_id,[representative_id])
            dist = distance_function(candidate_vec,representative_vec)
            x += 1
            if dist <= threshold
                specificity += 1
                if dist < closest_dist
                    closest_id = representative_id
                    closest_dist = dist
                end
            end
        end
        if isnothing(closest_id)
            unmapped_dict[candidate_id] = candidate_vec
        else
            push!(groups[closest_id],candidate_id)
            distances_dict[candidate_id] = closest_dist
            specificity_dict[candidate_id] = specificity
        end
    end
    return groups,unmapped_dict,distances_dict,specificity_dict,x
end

function optimization_process(;concat_representatives,merged_representatives,merged_groups,optimization_set,cache_dict,threshold,distance_function)
    n_comps = 0
    optimized_dict,n_comps_ = representative_optimization(
        representatives=merged_representatives,
        groups=merged_groups,
        optimization_set=optimization_set,
        cache_dict=cache_dict,
        distance_function=distance_function
    )
    n_comps += n_comps_

    optimized_groups,unmapped_dict,distances_dict,specificity_dict,n_comps_ = group_stratification_predefined_medoids(
        data_dict=concat_representatives,
        representative_dict=optimized_dict,
        cache_dict=cache_dict,
        threshold=threshold,
        distance_function=distance_function
    )
    n_comps += n_comps_

    shared_keys = intersect(keys(optimized_dict),keys(unmapped_dict))
    if !isempty(shared_keys)
        error("Overlapping object IDs in (optimized) representative_dict and unmapped_dict!")
    end
    merge!(optimized_dict,unmapped_dict)
    for id in keys(unmapped_dict)
        optimized_groups[id] = [id]
    end

    unmapped_set = Set(keys(unmapped_dict))

    return optimized_dict,optimized_groups,unmapped_set,distances_dict,specificity_dict,n_comps
end

function representative_optimization(;representatives,groups,optimization_set,cache_dict,distance_function)
    optimized_dict= Dict{String,Vector{Float32}}()
    n_comps = 0
    for (representative_id,group) in groups
        if representative_id in optimization_set
            if (length(group) == 1) # should never encounter if optimization_set is set properly -- remove after test
                error("Unexpected state: group ($representative_id) has a size of 1 -- unable to optimize!")
            end
            medoid_dict = get_medoid_dictionary_from_representatives(group,representative_dict,cache_dict)
            centroid_vec = compute_centroid(medoid_dict)
            medoid_id,medoid_vec,n_centroid_comps = determine_medoid(medoid_dict,centroid_vec,distance_function)
            n_comps += n_centroid_comps
            optimized_dict[medoid_id] = medoid_vec
        else
            optimized_dict[representative_id] = representatives[representative_id]
        end
    end
    return optimized_dict,n_comps
end
