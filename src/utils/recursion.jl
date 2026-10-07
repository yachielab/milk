module Recursion

    using Random
    using Logging
    using StatsBase

    using ..PairwiseComparisons: map_object_to_representative,map_object_to_representative_precomputed,
                                 determine_percentile_threshold,map_distance_function,compute_pairwise_distance_matrix
    using ..RepresentativeOptimization: optimize_representatives
    using ..GroupAggregation: group_aggregation
    using ..GroupStratification: stratification_predefined_medoids_precomputed_distances,stratification_precomputed_distances,
                                partitioned_group_stratification,stratification_process_direct_execution
    using ..FileHandling: partition_and_batch_input_files,clean_directory,load_milk_binaries,
                          write_group_results,attempt_to_load_cache

    export recursive_processing_framework,recursive_processing_direct_execution


    function recursive_processing_framework(;input_path,label,cache_path,invariant_args)

        partition_dir,files,batches = partition_and_batch_input_files(input_path,label,invariant_args)
        @info "\tInput file partitioned into $(length(files)) files ($(length(batches)) batches)"
        flush(stdout)

        @info "\tDetermining threshold: randomly sampling files..." 
        flush(stdout)
        if length(files) > 1
            if length(files) > invariant_args["batch-size"]
                Random.seed!(invariant_args["seed"])
                sampled_files = sample(files,invariant_args["batch-size"])
            else
                sampled_files = files
            end
            start_time = time()
            threshold = determine_percentile_threshold(sampled_files,invariant_args)
            elapsed_time = round((time()-start_time)/60,digits=2)
            @info "\tGlobally-determined percentile threshold: $(threshold) ($(length(sampled_files)) files; $elapsed_time min)"
            flush(stdout)
        else
            threshold = nothing
        end

        partitioned_group_stratification(
            batches=batches,
            files=files,
            threshold=threshold,
            cache_path=cache_path,
            label=label,
            partition_dir=partition_dir,
            invariant_args=invariant_args
        )

        representatives_path,groups_path = group_aggregation(
            label=label,
            threshold=threshold,
            cache_path=cache_path,
            partition_dir=partition_dir,
            invariant_args=invariant_args
        )
        exclusion_set = Set{String}([representatives_path])
        if !isnothing(cache_path) && isfile(cache_path)
            push!(exclusion_set,cache_path)
        end
        clean_directory(invariant_args["output-dir"],partition_dir,exclusion_set)
        return representatives_path,groups_path
    end


    function recursive_processing_direct_execution(;representatives_path,iteration,label,cache_path,invariant_args)
        """
        When population size is tractable (i.e., below partition threshold),
        carry out remaining recrusive iterations directly in a single function
        rather than the complete MILK framework, since the overhead becomes high
        for the latter.
        """
        @info "========== Starting encapsulated recursions..."

        distance_function = map_distance_function(invariant_args["metric"])

        data_dict,groupsize_dict = load_milk_binaries(representatives_path)
        n = length(data_dict)

        cache_dict = attempt_to_load_cache(cache_path)
        if isnothing(cache_dict)
            cache_dict = copy(data_dict)
        end

        i = iteration
        full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
        while n > invariant_args["sample-size"]
            @info "Iteration: $i ($n objects)"
            flush(stdout)
            data_dict,groupsize_dict = stratification_process_direct_execution(
                data_dict=data_dict,
                groupsize_dict=groupsize_dict,
                perc=invariant_args["percentile"],
                cache_dict=cache_dict,
                distance_function=distance_function,
                label=full_label,
                output_dir=invariant_args["output-dir"]
            )
            n = length(data_dict)
            @info "\t$n objects after recursive iteration."
            i += 1
            full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
        end
        return
    end

end

