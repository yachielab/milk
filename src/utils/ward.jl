module Ward

    using StatsBase

    export ward_groups,ward_percentile_threshold,compute_ward_cost_matrix

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

    function ward_groups(C,threshold,groupsizes)
        n = size(C,1)
        if n < 2
            return [[1]]
        end

        # tolerance: sqrt/square rounding of Float32 can land just below the cost it came from, blocking the final merge
        max_cost = threshold^2*(1+1f-5)

        sizes = copy(groupsizes)
        groups_dict = Dict(i => [i] for i in 1:n)
        active_groups = Set(1:n)

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
                c_update = ward_cost_update(C[i,k],C[j,k],c_ij,n_i,n_j,sizes[k])
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

end