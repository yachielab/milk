module GroupStratification

    using Glob
    using Distributed

    using ..FileHandling: attempt_to_load_cache,load_milk_binaries,write_milk_binaries,
                          write_group_results,open_file_write,
                          metadata_path,binary_path,groupsize_path
    using ..PairwiseComparisons: map_distance_function,compute_pairwise_distance_matrix,compute_euclidean_distance
    using ..HPC: submit_stratification_batch_array_jobs
    using ..Ward: compute_ward_cost_matrix,ward_groups,ward_percentile_threshold

    export stratification_process_direct_execution,
           stratification_process_distributed_execution,
           partitioned_group_stratification,
           batch_group_stratification_hpc_mode,
           ward_stratification,
           medoid_stratification

    function partitioned_group_stratification(;batches,files,threshold,cache_path,label,partition_dir,invariant_args)
        b = length(batches)
        p = invariant_args["percentile"]
        distance_function = map_distance_function(invariant_args["metric"])
        cache_dict        = attempt_to_load_cache(cache_path)
        if (b == 1 || !invariant_args["hpc-mode"])
            pathlist = sort(files)
            start_time = time()
            info_list = pmap(path -> stratification_process_distributed_execution(path,threshold,p,distance_function,invariant_args["linkage"],cache_dict,partition_dir), pathlist)
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
                cache_path=cache_path,
                label=label,
                partition_dir=partition_dir,
                invariant_args=invariant_args
            )
        end
        return
    end

    function batch_group_stratification_hpc_mode(;input_dir,threshold,percentile,metric,linkage,cache_path,verbose,output_dir)

        @info "Starting distributed stratification process"
        @info "  $(nworkers()) workers"
        @info "  Input batch directory: $(input_dir)"
        @info "  Threshold: $(threshold)"
        @info "  Percentile: $(percentile)"
        @info "  Metric: $(metric)"
        @info "  Linkage: $(linkage)"
        @info "  Cache path: $(cache_path)"
        @info "  Output directory: $(output_dir)"
        flush(stdout)

        n_comparisons = 0
        distance_function = map_distance_function(metric)
        cache_dict = attempt_to_load_cache(cache_path)

        @info "Initiating the pmap distributed call"
        flush(stdout)
        pathlist = sort(glob("*.ids",input_dir))
        info_list = pmap(path -> stratification_process_distributed_execution(path,threshold,percentile,distance_function,linkage,cache_dict,output_dir), pathlist)
        if verbose
            for info in info_list
                @info info
            end
        end

    end

    function stratification_process_direct_execution(;data_dict,groupsize_dict,perc,cache_dict,distance_function,linkage,label,output_dir)
        """
        This is only called when performing direct execution of recursive iterations.
        In such cases, there is only a single file to work with, so no need to distribute computing.
        """
        start_time = time()

        n_comparisons = 0
        D,index_map,threshold,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function)
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
    
        groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
        write_group_results(
            path=groups_path,
            label=label,
            stage="direct_stratification_process",
            cache_label="yes",
            n_input_objects=N,
            n_groups=G,
            groups=optimized_groups,
            direct_groupsize_dict=direct_groupsize_dict,
            total_groupsize_dict=total_groupsize_dict,
            distances_dict=distances_dict,
            specificity_dict=specificity_dict,
            threshold=threshold,
            n_comparisons=n_comparisons
        )
        return optimized_dict,total_groupsize_dict
    end

    function stratification_process_distributed_execution(path,threshold,perc,distance_function,linkage,cache_dict,output_dir)

        label = replace(basename(path),r"\.ids$" => "")

        start_time = time()

        n_comparisons = 0
        data_dict,groupsize_dict = load_milk_binaries(path)
        if length(data_dict) < 2
            representatives_path = joinpath(output_dir,"$(label).representatives.ids")
            groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
            if length(data_dict) == 0
                write_milk_binaries(Dict{String,Vector{Float32}}(),Dict{String,Int}(),representatives_path)
                open_file_write(groups_path) do file end
                touch(metadata_path(groups_path))

            else
                representative_id = collect(keys(data_dict))[1]
                write_milk_binaries(data_dict,groupsize_dict,representatives_path)
                write_group_results(
                    path=groups_path,
                    label=label,
                    stage="batch_stratification_process",
                    cache_label="no",
                    n_input_objects=length(data_dict),
                    n_groups=1,
                    groups=Dict(representative_id => [representative_id]),
                    direct_groupsize_dict=Dict(representative_id => 1),
                    total_groupsize_dict=Dict(representative_id => groupsize_dict[representative_id]),
                    distances_dict=Dict(representative_id => Float32[]),
                    specificity_dict=Dict(representative_id => Int32[]),
                    threshold=threshold,
                    n_comparisons=n_comparisons
                )
                rm(path)
                rm(binary_path(path))
                rm(groupsize_path(path))
            end
            info = "\t[$label] 1 group; $(length(data_dict) ) objects ($(length(data_dict) ) total); threshold: $threshold (predefined); 0 comparisons ($(nworkers()) CPU(s)); runtime: 0 seconds"
            touch(joinpath(output_dir,"$(label).flag"))
            return info
        else
            D,index_map,_,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function,false)
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
            cache_label = isnothing(cache_dict) ? "no" : "yes"
            elapsed_time = round((end_time-start_time)/60,digits=2)

            info = "\t[$label] $G groups; $n objects ($N total); threshold: $t (predefined); $n_comparisons comparisons ($(nworkers()) CPU(s)); runtime: $(elapsed_time) min"

            representatives_path = joinpath(output_dir,"$(label).representatives.ids")
            groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
            write_milk_binaries(optimized_dict,total_groupsize_dict,representatives_path)
            write_group_results(
                path=groups_path,
                label=label,
                stage="batch_stratification_process",
                cache_label=cache_label,
                n_input_objects=N,
                n_groups=G,
                groups=optimized_groups,
                direct_groupsize_dict=direct_groupsize_dict,
                total_groupsize_dict=total_groupsize_dict,
                distances_dict=distances_dict,
                specificity_dict=specificity_dict,
                threshold=threshold,
                n_comparisons=n_comparisons
            )
            rm(path)
            rm(binary_path(path))
            rm(groupsize_path(path))
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
