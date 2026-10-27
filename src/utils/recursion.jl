module Recursion

    using Random
    using Logging
    using StatsBase

    using ..PairwiseComparisons: determine_percentile_threshold,map_distance_function,compute_pairwise_distance_matrix
    using ..GroupAggregation: group_aggregation
    using ..GroupStratification: partitioned_group_stratification,stratification_process_direct_execution
    using ..FileHandling: partition_and_batch_input_files,clean_directory,load_milk_binaries,
                          write_group_results,load_milk_cache,load_group_spreads

    export recursive_processing_framework,
           recursive_processing_direct_execution


    function recursive_processing_framework(;input_path,label,invariant_args)

        partition_dir,files,batches = partition_and_batch_input_files(input_path,label,invariant_args)
        @info "\tInput file partitioned into $(length(files)) files ($(length(batches)) batches)"
        flush(stdout)

        @info "\tDetermining threshold: randomly sampling files..." 
        flush(stdout)
        if length(files) > 1
            if length(files) > invariant_args["batch-size"]
                Random.seed!(invariant_args["seed"])
                sampled_files = sample(files,invariant_args["batch-size"];replace=false)
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
            label=label,
            partition_dir=partition_dir,
            invariant_args=invariant_args
        )

        representatives_path,groups_path = group_aggregation(
            label=label,
            threshold=threshold,
            input_path=input_path,
            partition_dir=partition_dir,
            invariant_args=invariant_args
        )
        exclusion_set = Set{String}([representatives_path])
        clean_directory(invariant_args["output-dir"],partition_dir,exclusion_set)
        return representatives_path,groups_path
    end


    function recursive_processing_direct_execution(;representatives_path,iteration,label,invariant_args)
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

        cache_dict = invariant_args["cache-size"] >= 1 ? load_milk_cache(representatives_path,data_dict,groupsize_dict) : nothing
        spread_dict = load_group_spreads(representatives_path)
        i = iteration
        full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
        while n > invariant_args["sample-size"]
            @info "Iteration: $i ($n objects)"
            flush(stdout)
            n_previous = n
            data_dict,groupsize_dict,cache_dict,spread_dict = stratification_process_direct_execution(
                data_dict=data_dict,
                groupsize_dict=groupsize_dict,
                spread_dict=spread_dict,
                perc=invariant_args["percentile"],
                cache_dict=cache_dict,
                cache_size=invariant_args["cache-size"],
                distance_function=distance_function,
                linkage=invariant_args["linkage"],
                label=full_label,
                output_dir=invariant_args["output-dir"]
            )
            n = length(data_dict)

            if n == n_previous
                error("No groups merged at iteration $i ($n objects)! Stopping to avoid an infinity loop.")
            end

            @info "\t$n objects after recursive iteration."
            i += 1
            full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
        end
        return
    end

end

