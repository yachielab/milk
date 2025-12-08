
# include("pairwise_comparisons.jl")

function map_object(vec,representative_dict,threshold,distance_function)
    closest_id = nothing
    closest_dist = Inf
    specificity = 0
    n_comparisons = 0
    for (representative_id,representative_vec) in representative_dict
        dist = distance_function(vec,representative_vec)
        n_comparisons += 1
        if dist <= threshold
            specificity += 1
            if dist < closest_dist
                closest_id = representative_id
                closest_dist = dist
            end
        end
    end
    return closest_id,closest_dist,n_comparisons,specificity
end


function stratification(data_dict,threshold,distance_function)
    representative_dict = Dict{String,Vector{Float32}}()
    groups = Dict{String,Vector{String}}()
    x = 0
    for (candidate_id,candidate_vec) in data_dict
        closest_id,closest_dict,n_comps,_ = map_object(candidate_vec,representative_dict,threshold,distance_function)
        x += n_comps
        if isnothing(closest_id)
            representative_dict[candidate_id] = candidate_vec
            groups[candidate_id] = [candidate_id]
        else
            push!(groups[closest_id],candidate_id)
        end
    end
    return representative_dict,groups,x
end

function stratification_predefined_medoids(;data_dict,representative_dict,cache_dict,threshold,distance_function)
    groups = Dict( id => [id] for id in keys(representative_dict) )
    distances_dict = Dict(id => Vector{Float16}() for id in keys(representative_dict))
    specificity_dict = Dict(id => Vector{Int32}() for id in keys(representative_dict))
    x = 0
    for (candidate_id,candidate_vec) in data_dict
        if haskey(representative_dict,candidate_id) continue end
        closest_id,closest_dist,n_comps,specificity = map_object(candidate_vec,representative_dict,threshold,distance_function)
        x += n_comps
        if isnothing(closest_id)
            representative_dict[candidate_id] = candidate_vec
            groups[candidate_id] = [candidate_id]
            distances_dict[candidate_id] = Vector{Float16}()
            specificity_dict[candidate_id] = Vector{Int32}()
        else
            push!(groups[closest_id],candidate_id)
            push!(distances_dict[closest_id],closest_dist)
            push!(specificity_dict[closest_id],specificity)
        end
    end
    if !isnothing(cache_dict)
        for (cache_id,cache_vec) in cache_dict
            if haskey(data_dict,cache_id) continue end
            closest_id,closest_dist,n_comps,specificity = map_object(cache_vec,representative_dict,threshold,distance_function)
            x += n_comps
            if !isnothing(closest_id)
                push!(distances_dict[closest_id],closest_dist)
                push!(specificity_dict[closest_id],specificity)
            end
        end
    end
    return representative_dict,groups,distances_dict,specificity_dict,x
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
