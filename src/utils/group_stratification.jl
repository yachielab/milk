module GroupStratification

    using Glob
    using Distributed

    using ..FileHandling: load_input_array_as_dictionary,load_groups_as_dictionary,attempt_to_load_cache,
                          attempt_to_load_previous_groups,write_dictionary_as_csv,write_group_results
    using ..PairwiseComparisons: map_distance_function,map_object_to_representative_precomputed,
                                 map_object_to_representative,compute_pairwise_distance_matrix
    using ..RepresentativeOptimization: optimize_representatives

    export stratification,
           stratification_predefined_medoids,
           stratification_precomputed_distances,
           stratification_predefined_medoids_precomputed_distances,
           stratification_process_direct_execution,
           stratification_process_distributed_execution,
           compile_previous_groupings,
           partitioned_group_stratification

    function stratification_precomputed_distances(data_dict,D,index_map,threshold)
        representative_dict = Dict{String,Vector{Float32}}()
        groups = Dict{String,Vector{String}}()
        n_comparisons = 0
        for (candidate_id,candidate_vec) in data_dict
            i = index_map[candidate_id]
            closest_id,_,_ = map_object_to_representative_precomputed(i,representative_dict,D,index_map,threshold)
            if isnothing(closest_id)
                representative_dict[candidate_id] = candidate_vec
                groups[candidate_id] = [candidate_id]
            else
                push!(groups[closest_id],candidate_id)
            end
        end
        return representative_dict,groups,n_comparisons
    end

    function stratification(data_dict,threshold,distance_function)
        representative_dict = Dict{String,Vector{Float32}}()
        groups = Dict{String,Vector{String}}()
        n_comparisons = 0
        for (candidate_id,candidate_vec) in data_dict
            closest_id,_,n_comps,_ = map_object_to_representative(candidate_vec,representative_dict,threshold,distance_function)
            n_comparisons += n_comps
            if isnothing(closest_id)
                representative_dict[candidate_id] = candidate_vec
                groups[candidate_id] = [candidate_id]
            else
                push!(groups[closest_id],candidate_id)
            end
        end
        return representative_dict,groups,n_comparisons
    end

    function stratification_predefined_medoids_precomputed_distances(;data_dict,representative_dict,D,index_map,cache_dict,threshold,distance_function)
        groups = Dict( id => [id] for id in keys(representative_dict) )
        distances_dict = Dict(id => Vector{Float32}() for id in keys(representative_dict))
        specificity_dict = Dict(id => Vector{Int32}() for id in keys(representative_dict))
        for (candidate_id,candidate_vec) in data_dict
            if haskey(representative_dict,candidate_id) continue end
            i = index_map[candidate_id]
            closest_id,closest_dist,specificity = map_object_to_representative_precomputed(i,representative_dict,D,index_map,threshold)
            if isnothing(closest_id)
                """
                New groups can still be formed if candidates are not similar to any of the optimized medoids, where
                similarity is determined by the threshold.
                """
                representative_dict[candidate_id] = candidate_vec
                groups[candidate_id] = [candidate_id]
                distances_dict[candidate_id] = Vector{Float32}()
                specificity_dict[candidate_id] = Vector{Int32}()
            else
                push!(groups[closest_id],candidate_id)
                push!(distances_dict[closest_id],closest_dist)
                push!(specificity_dict[closest_id],specificity)
            end
        end
        x = 0
        if !isnothing(cache_dict)
            for (cache_id,cache_vec) in cache_dict
                if haskey(index_map,cache_id) continue end
                closest_id,closest_dist,n_comps,specificity = map_object_to_representative(cache_vec,representative_dict,threshold,distance_function)
                x += n_comps
                if !isnothing(closest_id)
                    push!(distances_dict[closest_id],closest_dist)
                    push!(specificity_dict[closest_id],specificity)
                end
            end
        end
        return representative_dict,groups,distances_dict,specificity_dict,x
    end

    function stratification_predefined_medoids(;data_dict,representative_dict,cache_dict,threshold,distance_function)
        groups = Dict( id => [id] for id in keys(representative_dict) )
        distances_dict = Dict(id => Vector{Float32}() for id in keys(representative_dict))
        specificity_dict = Dict(id => Vector{Int32}() for id in keys(representative_dict))
        x = 0
        for (candidate_id,candidate_vec) in data_dict
            if haskey(representative_dict,candidate_id) continue end
            closest_id,closest_dist,n_comps,specificity = map_object_to_representative(candidate_vec,representative_dict,threshold,distance_function)
            x += n_comps
            if isnothing(closest_id)
                representative_dict[candidate_id] = candidate_vec
                groups[candidate_id] = [candidate_id]
                distances_dict[candidate_id] = Vector{Float32}()
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
                closest_id,closest_dist,n_comps,specificity = map_object_to_representative(cache_vec,representative_dict,threshold,distance_function)
                x += n_comps
                if !isnothing(closest_id)
                    push!(distances_dict[closest_id],closest_dist)
                    push!(specificity_dict[closest_id],specificity)
                end
            end
        end
        return representative_dict,groups,distances_dict,specificity_dict,x
    end

    function partitioned_group_stratification(;batches,files,threshold,cache_path,previous_groups_path,label,partition_dir,invariant_args)
        b = length(batches)
        p = invariant_args["percentile"]
        distance_function = map_distance_function(invariant_args["metric"])
        cache_dict        = attempt_to_load_cache(cache_path)
        previous_groups   = attempt_to_load_previous_groups(previous_groups_path)
        if b == 1
            pathlist = sort(glob("*.csv",batches[1]))
            info_list = pmap(path -> stratification_process_distributed_execution(path,threshold,p,distance_function,cache_dict,previous_groups,partition_dir), pathlist)
            GC.gc() 
            if invariant_args["verbose"]
                for info in info_list
                    @info info
                end
            end
        else
            @info "\tSubmitting $b batches as an array job..."
            submit_stratification_batch_array_jobs(
                batches=batches,
                files=files,
                threshold=threshold,
                cache_path=cache_path,
                previous_groups_path=previous_groups_path,
                label=label,
                partition_dir=partition_dir,
                invariant_args=invariant_args
            )
        end
        return
    end

    function batch_group_stratification_hpc(input_dir,threshold,percentile,metric,cache_path,
                                  previous_groups_path,verbose,output_dir)

        @info "Starting distributed stratification process"
        @info "  $(nworkers()) workers"
        @info "  Input batch directory: $(input_dir)"
        @info "  Threshold: $(threshold)"
        @info "  Percentile: $(percentile)"
        @info "  Metric: $(metric)"
        @info "  Cache path: $(cache_path)"
        @info "  Previous groups path: $(previous_groups_path)"
        @info "  Output directory: $(output_dir)"
        flush(stdout)

        n_comparisons = 0
        distance_function = map_distance_function(metric)
        cache_dict = attempt_to_load_cache(cache_path)
        previous_groups = attempt_to_load_previous_groups(previous_groups_path)

        @info "Initiating the pmap distributed call"
        flush(stdout)
        pathlist = sort(glob("*.csv",input_dir))
        info_list = pmap(path -> stratification_process_distributed_execution(path,threshold,percentile,distance_function,cache_dict,previous_groups,output_dir), pathlist)

        if verbose
            for info in info_list
                @info info
            end
        end

    end

    function compile_previous_groupings(groups,previous_groups)
        compiled_groups = Dict{String,Vector{String}}()
        for (representative_id,group) in groups
            """
            This if-branch is checking for case when the previous and
            current groupings are exact same. No optimization of the medoid is required.
            """
            if (haskey(previous_groups,representative_id)) && (length(group) == 1)
                if group[1] != representative_id
                    error("Unexpected state: group[1] ($(group[1])) is not the representative_id ($representative_id)")
                end
                compiled_groups[representative_id] = previous_groups[representative_id]
            else
                compiled_groups[representative_id] = Vector{String}()
                for id in group
                    if !haskey(previous_groups, id)
                        error("Missing previous group for id: $id")
                    end
                    append!(compiled_groups[representative_id],previous_groups[id])
                end
            end
        end
        return compiled_groups
    end

    function stratification_process_direct_execution(;data_dict,perc,cache_dict,previous_groups,distance_function,label,output_dir)
        """
        This is only called when performing direct execution of recursive iterations.
        In such cases, there is only a single file to work with, so no need to distribute computing.
        """
        start_time = time()

        n_comparisons = 0
        D,index_map,threshold,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function)
        n_comparisons += x

        _,groups = stratification_precomputed_distances(data_dict,D,index_map,threshold)
        optimized_dict,optimization_set,x = optimize_representatives(data_dict,groups,cache_dict,distance_function)
        n_comparisons += x
        optimized_dict,optimized_groups,distances_dict,specificity_dict,x = stratification_predefined_medoids_precomputed_distances(
            data_dict=data_dict,
            representative_dict=optimized_dict,
            D=D,
            index_map=index_map,
            cache_dict=cache_dict,
            threshold=threshold,
            distance_function=distance_function
        )
        n_comparisons += x
        direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )

        previous_groups_label = "no"
        if !isnothing(previous_groups)
            previous_groups_label = "yes"
            optimized_groups = compile_previous_groupings(optimized_groups,previous_groups)
        end

        end_time = time()
        n = length(data_dict)
        G = length(optimized_dict)
        g = length(optimization_set)
        t = round(threshold,digits=4)
        N = sum(length(group) for group in values(optimized_groups))
        elapsed_time = round(end_time-start_time,digits=2)

        @info "\t[$label] $G groups ($g groups optimized); $n objects ($N total); threshold: $t (computed); $n_comparisons comparisons (1 CPU(s)); cache: yes; previous_groups: $previous_groups_label); runtime: $(elapsed_time) seconds"

        groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
        write_group_results(
            path=groups_path,
            label=label,
            stage="direct_stratification_process",
            cache_label="yes",
            compiled_label="yes",
            n_input_objects=N,
            n_groups=G,
            groups=optimized_groups,
            optimization_set=optimization_set,
            direct_groupsize_dict=direct_groupsize_dict,
            distances_dict=distances_dict,
            specificity_dict=specificity_dict,
            threshold=threshold,
            n_comparisons=n_comparisons
        )
        return optimized_dict,cache_dict,optimized_groups
    end

    function stratification_process_distributed_execution(path,threshold,perc,distance_function,cache_dict,previous_groups,output_dir)

        label = replace(basename(path),".csv" => "")

        start_time = time()

        n_comparisons = 0
        data_dict = load_input_array_as_dictionary(path)
        D,index_map,threshold_,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function)
        n_comparisons += x

        if isnothing(threshold)
            threshold_label = "computed"
            threshold = threshold_
        else
            threshold_label = "predefined"
        end

        _,groups = stratification_precomputed_distances(data_dict,D,index_map,threshold)
        optimized_dict,optimization_set,x = optimize_representatives(data_dict,groups,cache_dict,distance_function)
        n_comparisons += x
        optimized_dict,optimized_groups,distances_dict,specificity_dict,x = stratification_predefined_medoids_precomputed_distances(
            data_dict=data_dict,
            representative_dict=optimized_dict,
            D=D,
            index_map=index_map,
            cache_dict=cache_dict,
            threshold=threshold,
            distance_function=distance_function
        )
        n_comparisons += x
        direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )

        previous_groups_label = "no"
        if !isnothing(previous_groups)
            previous_groups_label = "yes"
            optimized_groups = compile_previous_groupings(optimized_groups,previous_groups)
        end

        end_time = time()

        n = length(data_dict)
        G = length(optimized_dict)
        g = length(optimization_set)
        t = round(threshold,digits=4)
        N = sum(length(group) for group in values(optimized_groups))
        cache_label = isnothing(cache_dict) ? "no" : "yes"
        elapsed_time = round(end_time-start_time,digits=2)

        info = "\t[$label] $G groups ($g groups optimized); $n objects ($N total); threshold: $t ($threshold_label); $n_comparisons comparisons ($(nworkers())) CPU(s); cache: $cache_label; previous_groups: $previous_groups_label); runtime: $(elapsed_time) seconds"

        representatives_path = joinpath(output_dir,"$(label).representatives.csv")
        groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
        write_dictionary_as_csv(optimized_dict,representatives_path)
        write_group_results(
            path=groups_path,
            label=label,
            stage="batch_stratification_process",
            cache_label=cache_label,
            compiled_label=previous_groups_label,
            n_input_objects=N,
            n_groups=G,
            groups=optimized_groups,
            optimization_set=optimization_set,
            direct_groupsize_dict=direct_groupsize_dict,
            distances_dict=distances_dict,
            specificity_dict=specificity_dict,
            threshold=threshold,
            n_comparisons=n_comparisons
        )
        rm(path)
        return info
    end


end
