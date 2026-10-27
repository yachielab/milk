module GroupStratification

    using Glob
    using Distributed

    using ..FileHandling: load_milk_binaries,write_milk_binaries,
                          write_group_results,open_file_write,
                          metadata_path,binary_path,group_meta_path,
                          GroupCache,cache_meta_path,cache_binary_path,count_rows,load_milk_cache
    using ..PairwiseComparisons: map_distance_function,compute_pairwise_distance_matrix,compute_euclidean_distance
    using ..HPC: submit_stratification_batch_array_jobs
    using ..Ward: compute_ward_cost_matrix,ward_groups,ward_percentile_threshold

    export stratification_process_direct_execution,
           stratification_process_distributed_execution,
           partitioned_group_stratification,
           batch_group_stratification_hpc_mode,
           ward_stratification,
           medoid_stratification,
           update_cache

    function partitioned_group_stratification(;batches,files,threshold,label,partition_dir,invariant_args)
        b = length(batches)
        p = invariant_args["percentile"]
        distance_function = map_distance_function(invariant_args["metric"])
        if (b == 1 || !invariant_args["hpc-mode"])
            pathlist = sort(files)
            start_time = time()
            info_list = pmap(path -> stratification_process_distributed_execution(path,threshold,p,distance_function,invariant_args["linkage"],invariant_args["cache-size"],partition_dir), pathlist)
            GC.gc()
            elapsed_time = round((time()-start_time)/60,digits=2)
            if invariant_args["verbose"]
                @info "\tPartitioned group stratification complete! $elapsed_time min to proceses $(length(files)) files"
            end
        else
            @info "\tSubmitting $b batches as an array job..."
            submit_stratification_batch_array_jobs(
                batches=batches,
                files=files,
                threshold=threshold,
                label=label,
                partition_dir=partition_dir,
                invariant_args=invariant_args
            )
        end
        return
    end

    function batch_group_stratification_hpc_mode(;input_dir,threshold,percentile,metric,linkage,cache_size,verbose,output_dir)

        @info "Starting distributed stratification process"
        @info "  $(nworkers()) workers"
        @info "  Input batch directory: $(input_dir)"
        @info "  Threshold: $(threshold)"
        @info "  Percentile: $(percentile)"
        @info "  Metric: $(metric)"
        @info "  Linkage: $(linkage)"
        @info "  Cache size: $(cache_size)"
        @info "  Output directory: $(output_dir)"
        flush(stdout)

        distance_function = map_distance_function(metric)

        @info "Initiating the pmap distributed call"
        flush(stdout)
        pathlist = sort(glob("*.ids",input_dir))
        info_list = pmap(path -> stratification_process_distributed_execution(path,threshold,percentile,distance_function,linkage,cache_size,output_dir), pathlist)
        if verbose
            for info in info_list
                @info info
            end
        end

    end

    function stratification_process_direct_execution(;data_dict,groupsize_dict,perc,cache_dict,cache_size,distance_function,linkage,label,output_dir)
        """
        This is only called when performing direct execution of recursive iterations.
        In such cases, there is only a single file to work with, so no need to distribute computing.
        """
        start_time = time()

        n_comparisons = 0
        D,index_map,threshold,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function;cache_dict=(cache_size >= 1 ? cache_dict : nothing))
        n_comparisons += x
        if startswith(linkage,"ward")
            balanced = linkage == "ward_balanced"
            optimized_dict,optimized_groups,distances_dict,threshold = ward_stratification(data_dict,groupsize_dict,D,index_map,nothing,perc,balanced)
            specificity_dict = Dict(id => Int32[] for id in keys(optimized_groups))
        else
            optimized_dict,optimized_groups,distances_dict,specificity_dict = medoid_stratification(data_dict,groupsize_dict,D,index_map,threshold)
        end
        direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )
        total_groupsize_dict = Dict(id => sum(groupsize_dict[m] for m in group) for (id,group) in optimized_groups)


        end_time = time()
        n = length(data_dict)
        G = length(optimized_dict)
        t = round(threshold,digits=4)
        N = sum(values(total_groupsize_dict))
        elapsed_time = round((end_time-start_time)/60,digits=2)

        @info "\t[$label] $G groups; $n objects ($N total); threshold: $t (computed); $n_comparisons comparisons (1 CPU(s)); runtime: $(elapsed_time) min"
        
        cache_dict = cache_size >= 1 ? update_cache(optimized_dict,optimized_groups,cache_dict,cache_size,distance_function) : nothing
        groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
        write_group_results(
            path=groups_path,
            label=label,
            stage="direct_stratification_process",
            n_cached=isnothing(cache_dict) ? G : sum(length(cache.cell_ids) for cache in values(cache_dict)),
            n_input_objects=N,
            n_groups=G,
            groups=optimized_groups,
            direct_groupsize_dict=direct_groupsize_dict,
            total_groupsize_dict=total_groupsize_dict,
            distances_dict=distances_dict,
            specificity_dict=specificity_dict,
            threshold=threshold,
            n_comparisons=n_comparisons,
            cache_dict=cache_dict
        )
        return optimized_dict,total_groupsize_dict,cache_dict
    end

    """
    Cache of one new group: its representative plus up to `cache_size` extra cells, picked
    farthest-first from the cached cells of the groups it merged. Each cached cell stands for the
    pooled cells nearest to it (subgroup_size = their summed sizes, subgroup_radius = the largest
    distance to them plus their own radius).
    """
    function _update_cache_group(representative_id,representative_vec,member_caches,cache_size,distance_function)
        pooled_ids = reduce(vcat,(cache.cell_ids for cache in member_caches))
        pooled_vecs = reduce(vcat,(cache.vecs for cache in member_caches))
        pooled_sizes = reduce(vcat,(cache.subgroup_sizes for cache in member_caches))
        pooled_radii = reduce(vcat,(cache.subgroup_radii for cache in member_caches))
        n = length(pooled_ids)

        cache_ids = [representative_id]
        cache_vecs = [representative_vec]

        # each pooled cell -> its nearest cached cell
        nearest_dist = [distance_function(vec,representative_vec) for vec in pooled_vecs]
        nearest_cached = ones(Int,n) # initially, they are all closest to representative
        nearest_dist[pooled_ids .== representative_id] .= 0 # representative's own row is already cached

        # pos=1 is representative; plus one since "in addition to" representative
        for k in 2:cache_size+1
            j = argmax(nearest_dist) # farthest pooled cell from the cache
            nearest_dist[j] == 0 && break # everything already cached
            push!(cache_ids,pooled_ids[j])
            push!(cache_vecs,pooled_vecs[j])
            for (i,vec) in enumerate(pooled_vecs)
                updated_dist = distance_function(vec,pooled_vecs[j])
                if updated_dist < nearest_dist[i]
                    nearest_dist[i] = updated_dist
                    nearest_cached[i] = k
                end
            end
            nearest_dist[j] = 0 # exact: a cell's computed distance to itself can be ~1e-7 (cosine rounding)
            nearest_cached[j] = k
        end

        subgroup_sizes = zeros(Int,length(cache_ids))
        subgroup_radii = zeros(Float64,length(cache_ids))
        for i in 1:n
            subgroup_sizes[nearest_cached[i]] += pooled_sizes[i]
            subgroup_radii[nearest_cached[i]] = max(subgroup_radii[nearest_cached[i]],nearest_dist[i]+pooled_radii[i])
        end

        return GroupCache(cache_ids,cache_vecs,subgroup_sizes,subgroup_radii)
    end

    function update_cache(representative_dict,groups_dict,cache_dict,cache_size,distance_function)
        # note that members sorted so pooling order (and tie-breaking) doesn't depend on Dict order
        updated_cache_dict = Dict{String,GroupCache}()
        for (id,vec) in representative_dict
            member_caches = [cache_dict[m] for m in sort(groups_dict[id])]
            updated_cache_dict[id] =  _update_cache_group(id,vec,member_caches,cache_size,distance_function)
        end
        return updated_cache_dict
    end

    function stratification_process_distributed_execution(path,threshold,perc,distance_function,linkage,cache_size,output_dir)

        label = replace(basename(path),r"\.ids$" => "")

        start_time = time()

        n_comparisons = 0
        data_dict,groupsize_dict = load_milk_binaries(path)
        cache_dict = cache_size >= 1 ? load_milk_cache(path,data_dict,groupsize_dict) : nothing
        if length(data_dict) < 2
            representatives_path = joinpath(output_dir,"$(label).representatives.ids")
            groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
            if length(data_dict) == 0
                write_milk_binaries(Dict{String,Vector{Float32}}(),Dict{String,Int}(),representatives_path;cache_dict=cache_dict)
                open_file_write(groups_path) do file end
                touch(metadata_path(groups_path))
            else
                representative_id = collect(keys(data_dict))[1]
                write_milk_binaries(data_dict,groupsize_dict,representatives_path;cache_dict=cache_dict)
                write_group_results(
                    path=groups_path,
                    label=label,
                    stage="batch_stratification_process",
                    n_cached=count_rows(representatives_path),
                    n_input_objects=length(data_dict),
                    n_groups=1,
                    groups=Dict(representative_id => [representative_id]),
                    direct_groupsize_dict=Dict(representative_id => 1),
                    total_groupsize_dict=Dict(representative_id => groupsize_dict[representative_id]),
                    distances_dict=Dict(representative_id => Float32[]),
                    specificity_dict=Dict(representative_id => Int32[]),
                    threshold=threshold,
                    n_comparisons=n_comparisons,
                    cache_dict=cache_dict
                )
                rm(path)
                rm(binary_path(path))
                rm(group_meta_path(path))
                rm(cache_meta_path(path),force=true)
                rm(cache_binary_path(path),force=true)
            end
            info = "\t[$label] 1 group; $(length(data_dict) ) objects ($(length(data_dict) ) total); threshold: $threshold (predefined); 0 comparisons ($(nworkers()) CPU(s)); runtime: 0 seconds"
            touch(joinpath(output_dir,"$(label).flag"))
            return info
        else
            D,index_map,_,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function,false;cache_dict=cache_dict)
            n_comparisons += x
            if startswith(linkage,"ward")
                balanced = linkage == "ward_balanced"
                optimized_dict,optimized_groups,distances_dict,threshold = ward_stratification(data_dict,groupsize_dict,D,index_map,threshold,perc,balanced)
                specificity_dict = Dict(id => Int32[] for id in keys(optimized_groups))
            else
                optimized_dict,optimized_groups,distances_dict,specificity_dict = medoid_stratification(data_dict,groupsize_dict,D,index_map,threshold)
            end
            direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )
            total_groupsize_dict = Dict( id => sum(groupsize_dict[m] for m in group) for (id,group) in optimized_groups )
            end_time = time()

            n = length(data_dict)
            G = length(optimized_dict)
            t = round(threshold,digits=4)
            N = sum(values(total_groupsize_dict))
            elapsed_time = round((end_time-start_time)/60,digits=2)

            info = "\t[$label] $G groups; $n objects ($N total); threshold: $t (predefined); $n_comparisons comparisons ($(nworkers()) CPU(s)); runtime: $(elapsed_time) min"

            representatives_path = joinpath(output_dir,"$(label).representatives.ids")
            groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
            cache_dict = cache_size >= 1 ? update_cache(optimized_dict,optimized_groups,cache_dict,cache_size,distance_function) : nothing
            write_milk_binaries(optimized_dict,total_groupsize_dict,representatives_path;cache_dict=cache_dict)
            write_group_results(
                path=groups_path,
                label=label,
                stage="batch_stratification_process",
                n_cached=count_rows(representatives_path),
                n_input_objects=N,
                n_groups=G,
                groups=optimized_groups,
                direct_groupsize_dict=direct_groupsize_dict,
                total_groupsize_dict=total_groupsize_dict,
                distances_dict=distances_dict,
                specificity_dict=specificity_dict,
                threshold=threshold,
                n_comparisons=n_comparisons,
                cache_dict=cache_dict
            )
            rm(path)
            rm(binary_path(path))
            rm(group_meta_path(path))
            rm(cache_meta_path(path),force=true)
            rm(cache_binary_path(path),force=true)
            touch(joinpath(output_dir,"$(label).flag"))
            return info
        end
    end

    function ward_stratification(data_dict,groupsize_dict,D,index_map,threshold,p,balanced)
        ids = Vector{String}(undef,length(index_map))
        for (id,i) in index_map
            ids[i] = id
        end

        sizes = [groupsize_dict[id] for id in ids]
        cost_sizes = balanced ? ones(Int,length(sizes)) : copy(sizes)
        C = compute_ward_cost_matrix(D,cost_sizes)
        if isnothing(threshold)
            threshold = ward_percentile_threshold(C,p)
        end
        groups = ward_groups(C,threshold,cost_sizes)

        representative_dict = Dict{String,Vector{Float32}}()
        groups_dict = Dict{String,Vector{String}}()
        distances_dict = Dict{String,Vector{Float32}}()
        for group in groups
            weights = sizes[group]
            centroid = sum(weights[x]*data_dict[ids[g]] for (x,g) in enumerate(group)) / sum(weights)
            dists = [compute_euclidean_distance(data_dict[ids[g]],centroid) for g in group]
            representative_id = ids[group[argmin(dists)]] # named after the member closest to the centroid

            representative_dict[representative_id] = centroid
            groups_dict[representative_id] = ids[group]
            distances_dict[representative_id] = dists
        end
        
        return representative_dict,groups_dict,distances_dict,threshold
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

        # soft-kernel density ordering of cells (most->least dense)
        density_vec = _soft_kernel_density(D,sizes,threshold)
        ordered_indices = sortperm(density_vec,rev=true,alg=MergeSort)

        # key: medoid index; value: vector of member indices
        medoid_dict = Dict{Int,Vector{Int}}()
        for iter in 1:max_iter

            updated_medoid_dict = Dict{Int,Vector{Int}}(medoid => Int[] for medoid in keys(medoid_dict))
            update_flag = false
            for i in ordered_indices
                if all(D[i,medoid]+radii[i] > threshold for medoid in keys(updated_medoid_dict))
                    updated_medoid_dict[i] = Int[]
                    update_flag = true
                end
            end
    
            for i in 1:n
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
        return medoid_dict
    end

    function medoid_stratification(data_dict,groupsize_dict,D,index_map,threshold;radius_dict=nothing)
        ids = Vector{String}(undef,length(index_map))
        for (id,i) in index_map
            ids[i] = id
        end

        sizes = Float64[groupsize_dict[id] for id in ids]
        radii = isnothing(radius_dict) ? zeros(length(ids)) : Float64[radius_dict[id] for id in ids]
        medoid_dict = threshold_bound_kmedoids(D,sizes,threshold,radii=radii)
        medoids = collect(keys(medoid_dict))

        representative_dict = Dict{String,Vector{Float32}}()
        groups_dict = Dict{String,Vector{String}}()
        distances_dict = Dict{String,Vector{Float32}}()
        specificity_dict = Dict{String,Vector{Int32}}()
        for (medoid,members) in medoid_dict
            medoid_id = ids[medoid]
            representative_dict[medoid_id] = data_dict[medoid_id]
            groups_dict[medoid_id] = ids[members]
            distances_dict[medoid_id] = Float32[D[i,medoid] for i in members if i != medoid]
            specificity_dict[medoid_id] = Int32[count(D[i,m]+radii[i] <= threshold for m in medoids) for i in members if i != medoid]
        end
        return representative_dict,groups_dict,distances_dict,specificity_dict
    end
end
