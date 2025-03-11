using Distances
using Logging

function compute_cosine_distance(vec_i,vec_j)
    return cosine_dist(vec_i,vec_j)
end

function compute_euclidean_distance(vec_i,vec_j)
    L = length(vec_i)
    d = 0
    @inbounds @simd for x in eachindex(vec_i,vec_j)
        d += (vec_i[x]-vec_j[x])^2
    end
    return sqrt(d)
end

function compute_manhattan_distance(vec_i,vec_j)
    L = length(vec_i)
    d = 0
    @inbounds @simd for x in eachindex(vec_i,vec_j)
        d += abs(vec_i[x]-vec_j[x])
    end
    return d
end

function map_distance_function(metric)
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

# Helper function
function pairwise_distance_thread(profiles,indices,distance_function)
    results = Vector{Tuple{UInt16,UInt16,Float32}}(undef,length(indices))
    for (x,(i,j)) in enumerate(indices)
        results[x] = (i,j,distance_function(profiles[i],profiles[j]))
    end
    return results
end

function pairwise_distances(profiles,threads,perc,distance_function)

    n = length(profiles)
    indices = Vector{Tuple{UInt16,UInt16}}([(i,j) for i in 1:n-1 for j in i+1:n])
    L = length(indices)

    # Multithreaded computation of pairwise distances
    partitions = Iterators.partition(indices,ceil(Int,L/threads))
    tasks = [@spawn pairwise_distance_thread(profiles,p,distance_function) for p in partitions]
    results = fetch.(tasks)

    D = zeros(Float32,n,n)
    distances = Vector{Float32}()
    for result in results
        for (i,j,d) in result
            D[i,j] = d
            D[j,i] = d
            push!(distances,d)
        end
    end

    n_comp = length(distances)
    threshold = percentile(distances,perc)

    # @info "\tThreshold: $threshold"

    return D,threshold,n_comp
end
