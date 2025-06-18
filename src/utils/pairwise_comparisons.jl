using Logging
using Distributed

@everywhere using StatsBase, Distances
@everywhere include("file_handling.jl")

@everywhere function compute_cosine_distance(vec_i,vec_j)
    return cosine_dist(vec_i,vec_j)
end

@everywhere function compute_euclidean_distance(vec_i,vec_j)
    L = length(vec_i)
    d = 0
    @inbounds @simd for x in eachindex(vec_i,vec_j)
        d += (vec_i[x]-vec_j[x])^2
    end
    return sqrt(d)
end

@everywhere function compute_manhattan_distance(vec_i,vec_j)
    L = length(vec_i)
    d = 0
    @inbounds @simd for x in eachindex(vec_i,vec_j)
        d += abs(vec_i[x]-vec_j[x])
    end
    return d
end

@everywhere function map_distance_function(metric)
    if metric == "cosine"
        distance_function = compute_cosine_distance
    elseif metric == "euclidean"
        distance_function = compute_euclidean_distance
    elseif metric == "manhattan"
        distance_function = compute_manhattan_distance
    else
        throw("Specified metric ($metric) not supported!")
    end
    return distance_function
end

@everywhere function compute_pairwise_distance_matrix(data_dict,perc,distance_function)
    """
    When working with lower-dimensional representations of gene expression space,
    distribution of distance matrix computation across files with single CPU distance matrix 
    computation is much faster than sequentially computing distance matrices with multiple
    threads.
    """
    ids = collect(keys(data_dict))
    index_map = Dict( ids[i] => i for i in eachindex(ids) )
    n = length(ids)
    D = zeros(Float16,n,n)
    indices = [ (ids[i],ids[j]) for i in 1:n-1 for j in i+1:n ]
    L = length(indices)
    distances = Vector{Float16}(undef,L)
    @inbounds @simd for x in eachindex(indices)
        id_i,id_j = indices[x]
        i = index_map[id_i]
        j = index_map[id_j]
        d = distance_function(data_dict[id_i],data_dict[id_j])
        D[i,j] = d
        D[j,i] = d
        distances[x] = d
    end
    threshold = percentile(distances,perc)
    x = length(distances)
    return D,index_map,threshold,x
end

@everywhere function determine_percentile_threshold_process(input_path,metric,perc)
    data_dict = load_input_array_as_dictionary(input_path)
    distance_function = map_distance_function(metric)
    D,index_map,threshold,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function)
    return threshold
end

function determine_percentile_threshold_distributed(pathlist,invariant_args)
    metric = invariant_args["metric"]
    perc = invariant_args["percentile"]
    thresholds = pmap(path -> determine_percentile_threshold_process(path,metric,perc),pathlist)
    @everywhere GC.gc() 
    return minimum(thresholds)
end
