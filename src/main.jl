function main()
    args = parse_arguments()

    if args["group-stratification-mode"]
        batch_group_stratification_hpc(
            input_dir=args["stratification-input-dir"],
            threshold=args["stratification-threshold"],
            percentile=args["stratification-percentile"],
            metric=args["stratification-metric"],
            cache_path=args["stratification-cache-path"],
            previous_groups_path=args["stratification-previous-groups-path"],
            verbose=args["stratification-verbose"],
            output_dir=args["stratification-output-dir"]
        )
    else
        logger = ConsoleLogger(stdout,Logging.Info)
        global_logger(logger)

        @info "MILK"
        flush(stdout)

        @info "Using $(nworkers()) worker(s) for distributed processing"

        log_args(args)
        invariant_args = instantiate_invariant_args(args)

        if isdir(args["output-dir"])
            if args["force-overwrite"]
                @warn "Output directory already exists! Overwriting."
                rm(args["output-dir"],recursive=true)
            else
                error("Output directory ($(args["output-dir"])) already exists! Exiting.")
            end
        end
        mkdir(invariant_args["output-dir"])

        if isnothing(args["label"])
            label = replace(basename(args["input-path"]),".csv" => "")
        else
            label = args["label"]
        end

        cache_path = nothing
        previous_groups_path = nothing
        i = 0 # Initial recursion

        full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"

        input_path = joinpath(invariant_args["output-dir"],"$(full_label).input.csv")
        absolute_input_path = isabspath(args["input-path"]) ? args["input-path"] : joinpath(pwd(),args["input-path"])
        if islink(input_path)
            @warn "Symlink already exists! ($input_path)"
            rm(input_path)
        end
        symlink(absolute_input_path,input_path)

        n = get_object_count(input_path)

        @info "Iteration: $i ($n objects)"
        flush(stdout)

        # Pre-emptively set cache for subsequent recursions (no caching for initial iteration)
        cache_path = attempt_to_cache_file(input_path,n,invariant_args)

        if n <= args["partition-size"]
            recursive_processing_direct_execution(
                representatives_path=input_path,
                iteration=i,
                label=label,
                cache_path=cache_path,
                previous_groups_path=previous_groups_path,
                invariant_args=invariant_args
            )
        else
            representatives_path,groups_path = recursive_processing_framework(
                input_path=input_path,
                label=full_label,
                cache_path=nothing,
                previous_groups_path=nothing,
                invariant_args=invariant_args
            )

            n = get_object_count(representatives_path)
            if isnothing(cache_path)
                cache_path = attempt_to_cache_file(representatives_path,n,invariant_args)
            end
            # previous_groups_path = groups_path
            if n <= args["compile-previous-threshold"]
                previous_groups_path = groups_path
            end

            while n > args["sample-size"]
                if n <= args["partition-size"]
                    recursive_processing_direct_execution(
                        representatives_path=representatives_path,
                        iteration=i,
                        label=label,
                        cache_path=cache_path,
                        previous_groups_path=previous_groups_path,
                        invariant_args=invariant_args
                    )
                    break
                end

                i += 1
                full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
                @info "Iteration: $i ($n objects)"
                flush(stdout)

                input_path = joinpath(invariant_args["output-dir"],"$(full_label).input.csv")
                symlink(representatives_path,input_path)

                representatives_path,groups_path = recursive_processing_framework(
                    input_path=input_path,
                    label=full_label,
                    cache_path=cache_path,
                    previous_groups_path=previous_groups_path,
                    invariant_args=invariant_args
                )
                # previous_groups_path = groups_path
                n = get_object_count(representatives_path)
                if isnothing(cache_path)
                    cache_path = attempt_to_cache_file(representatives_path,n,invariant_args)
                end
                if n <= args["compile-previous-threshold"]
                    previous_groups_path = groups_path
                end
                @info "\t$n objects after recursive iteration."
                flush(stdout)
            end
        end

        @info "\nCompleted recursive downsampling procedure."
        flush(stdout)

        if args["skip-reconstruction"]
            @info "Skipping cell hierarchy reconstruction."
        else
            @info "Reconstructing hierarchical graph..."
            hierarchical_reconstruction(args["output-dir"])
        end
        @info "\tDone!"
        flush(stdout)
    end
end
