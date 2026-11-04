module Linkage

    using StatsBase

    export ward_groups,
           ward_percentile_threshold,
           compute_ward_cost_matrix,
           compute_average_cost_matrix,
           average_cost_update,
           identify_connected_objects,
           threshold_bound_kmedoids

    """
    Objects that can group with at least one other: some j with C[i,j] + radii[i] <= threshold (i could join j),
    checked in both directions. Every other object can only ever be a singleton, so each linkage skips them:
        medoid: nothing is within its reach and it is within nobody's reach
        ward/average: reducible linkages, so merging other groups never lowers its costs
    With a low threshold this removes most objects of a partition.
    """
    function identify_connected_objects(C,threshold;radii=nothing)
        n = size(C,1)
        if isnothing(radii)
            radii = zeros(n)
        end
        connected_list = falses(n)
        @inbounds for j in 1:n
            for i in 1:n
                if isnan(C[i,j])
                    error("NaN distance between objects $(i) and $(j)")
                end
                if (i != j) && (C[i,j]+radii[i] <= threshold)
                    connected_list[i] = true
                    connected_list[j] = true
                end
            end
        end
        return connected_list
    end

    """
    Lance-Williams formula, Ward linkage: squared merge cost update so you don't have to recompute cost
        matrix every time a group merge occurs (https://en.wikipedia.org/wiki/Ward%27s_method)
        c refers to squared pairwise (merge) cost
        n refers to group size
    """
    function ward_cost_update(c_ik,c_jk,c_ij,n_i,n_j,n_k)
        N = n_i+n_j+n_k
        a_i = (n_i+n_k)/N
        a_j = (n_j+n_k)/N
        b = n_k/N
        return max(a_i*c_ik+a_j*c_jk-b*c_ij,0) # clamp for floating errors 
    end

    """
    Squared Ward merge cost (i.e., distance) from squared centroid distance d
    Factor of 2 makes the cost equal d when merging two singletons (matches scipy)
    """
    function squared_ward_cost(d,n_i,n_j)
        a = 2*(n_i*n_j)/(n_i+n_j)
        return a*d
    end

    function ward_percentile_threshold(C,p)
        n = size(C,1)
        return sqrt(percentile([C[i,j] for j in 2:n for i in 1:j-1],p))
    end

    function compute_ward_cost_matrix(D,groupsizes)
        n = size(D,1)
        C = zeros(Float32,n,n)
        for (i,j) in ((i,j) for j in 2:n for i in 1:j-1)
            c = squared_ward_cost(D[i,j]^2,groupsizes[i],groupsizes[j])
            C[i,j] = c
            C[j,i] = c
        end
        return C
    end


    """
    Average linkage on squared Euclidean distance (UPGMA on d^2):
        cost(i,j) = mean squared distance over all cell pairs between groups i and j
                  = squared centroid distance + spread_i + spread_j
    spread = mean squared distance of a group's cells to its centroid (0 for a single cell)
    """
    function compute_average_cost_matrix(D,spreads)
        n = size(D,1)
        C = zeros(Float32,n,n)
        for j in 2:n
            for i in 1:j-1
                c = D[i,j]^2 + spreads[i] + spreads[j]
                C[i,j] = c
                C[j,i] = c
            end
        end
        return C
    end

    """
    Lance-Williams update for average linkage: cost to a merged group is the size-weighted
    mean of the costs to its two parts (c_ij and n_k unused; same signature as ward_cost_update)
    """
    function average_cost_update(c_ik,c_jk,c_ij,n_i,n_j,n_k)
        return (n_i*c_ik + n_j*c_jk)/(n_i+n_j)
    end

    function ward_groups(C,threshold,groupsizes;update=ward_cost_update)
        n = size(C,1)
        if n < 2
            return [[1]]
        end

        # tolerance: sqrt/square rounding of Float32 can land just below the cost it came from, blocking the final merge
        max_cost = threshold^2*(1+1f-5)

        sizes = copy(groupsizes)
        groups_dict = Dict(i => [i] for i in 1:n)
        active_groups = Set(findall(identify_connected_objects(C,max_cost)))

        #=
        Nearest-neighbor chain: avoids searching the matrix for the global cheapest pair after every merge
            - chain of groups, each the cheapest partner of the one before: C(g_1,g_2) > C(g_2,g_3) > ...
            - push a partner only if its cost is below threshold
            - merge and pop g_end and g_penultimate when they are mutually cheapest
            - retire g_end if its cheapest cost > threshold (only possible when the chain is [g_end])
        =#
        chain = Int[] # stack-like behavior
        while length(active_groups) > 1

            if isempty(chain)
                push!(chain,first(active_groups))
            end

            i = chain[end]
            previous = length(chain) >= 2 ? chain[end-1] : nothing

            j = previous
            min_cost = isnothing(previous) ? Inf32 : C[i,j]

            for k in active_groups
                if k == i
                    continue
                end
                if C[k,i] < min_cost
                    j = k
                    min_cost = C[k,i]
                end
            end

            if min_cost > max_cost
                pop!(chain) # last element
                delete!(active_groups,i)
                continue
            elseif j != previous
                push!(chain,j)
                continue
            end

            # MERGE found: i and j are mutually cheapest partners
            pop!(chain)
            pop!(chain)

            n_i = sizes[i]
            n_j = sizes[j]
            N = n_i+n_j
            c_ij = C[i,j]
            for k in active_groups
                if k == i || k == j
                    continue
                end
                c_update = update(C[i,k],C[j,k],c_ij,n_i,n_j,sizes[k])
                C[i,k] = c_update
                C[k,i] = c_update
            end

            # merging j into i (arbitrary)
            sizes[i] = N
            append!(groups_dict[i],groups_dict[j])
            pop!(groups_dict,j)
            delete!(active_groups,j)
        end

        return collect(values(groups_dict))
    end

function _soft_kernel_density(D,sizes,threshold)
        n = size(D,1)
        bandwidth = threshold > 0 ? threshold : eps(Float64)
        density_vec = Vector{Float64}(undef,n)
        for j in 1:n
            kernels = exp.(-0.5 .* (view(D,:,j) ./ bandwidth).^2)
            density_vec[j] = sum(sizes .* kernels)
        end
        return density_vec
    end

    function _determine_weighted_medoid(D,sizes,members,density_vec)
        candidate_costs = Vector{Float64}(undef,length(members))
        for (k,candidate) in enumerate(members)
            candidate_costs[k] = sum(sizes[i]*D[i,candidate] for i in members)
        end
        min_cost = minimum(candidate_costs)
        tolerance = 1e-6*max(min_cost,1e-12)
        tied_candidates = members[candidate_costs .<= min_cost + tolerance]
        return tied_candidates[argmax(density_vec[tied_candidates])]
    end

    function threshold_bound_kmedoids(D,sizes,threshold;max_iter=50,radii=nothing)
        n = size(D,1)

        if isnothing(radii)
            radii = zeros(n)
        end

        start_time = time()

        # soft-kernel density ordering of cells (most->least dense)
        density_vec = _soft_kernel_density(D,sizes,threshold)
        ordered_indices = sortperm(density_vec,rev=true,alg=MergeSort)

        connected_list = identify_connected_objects(D,threshold;radii=radii)
        connected_indices = findall(connected_list)
        ordered_indices = [i for i in ordered_indices if connected_list[i]] # same density order
    
        # key: medoid index; value: vector of member indices
        medoid_dict = Dict{Int,Vector{Int}}()

        n_iterations = 0
        for iter in 1:max_iter
            n_iterations = iter

            updated_medoid_dict = Dict{Int,Vector{Int}}(medoid => Int[] for medoid in keys(medoid_dict))
            update_flag = false
            for i in ordered_indices
                if all(D[i,medoid] + radii[i] > threshold for medoid in keys(updated_medoid_dict))
                    updated_medoid_dict[i] = Int[]
                    update_flag = true
                end
            end
    
            for i in connected_indices
                if haskey(updated_medoid_dict,i)
                    push!(updated_medoid_dict[i],i)
                    continue
                end
                nearest_medoid = 0 # invalid index
                nearest_dist = Inf32
                for medoid in keys(updated_medoid_dict)
                    if D[i,medoid] + radii[i] <= threshold && D[i,medoid] < nearest_dist
                        nearest_medoid = medoid
                        nearest_dist = D[i,medoid]
                    end
                end
                nearest_medoid == 0 && error("object $i has no medoid within threshold (NaN distance?)")
                push!(updated_medoid_dict[nearest_medoid],i)
            end

            if (!update_flag && updated_medoid_dict == medoid_dict) || iter == max_iter
                medoid_dict = updated_medoid_dict
                break
            end

            medoid_dict = Dict{Int,Vector{Int}}()
            for (medoid,members) in updated_medoid_dict
                optimized_medoid = _determine_weighted_medoid(D,sizes,members,density_vec)
                medoid_dict[optimized_medoid] = members
            end
        end
        for (i,connected_flag) in enumerate(connected_list)
            if !connected_flag
                medoid_dict[i] = [i]
            end
        end
        elapsed_time = round(time()-start_time,digits=2)
        converged = n_iterations < max_iter
        @info "\tkmedoids: $(n-length(connected_indices)) isolated, $(length(connected_indices)) grouped over $(n_iterations) iterations ($(converged ? "converged" : "hit max_iter")); $(length(medoid_dict)) medoids, $(n) objects; $(elapsed_time)s"
        return medoid_dict
    end

end