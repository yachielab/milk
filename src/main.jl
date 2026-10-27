# function main()
function main()
    args = parse_arguments()
    validate_args(args)

    if args["group-stratification-mode"]
        batch_group_stratification_hpc_mode(
            input_dir=args["stratification-input-dir"],
            threshold=args["stratification-threshold"],
            percentile=args["stratification-percentile"],
            metric=args["stratification-metric"],
            linkage=args["stratification-linkage"],
            cache_size=args["stratification-cache-size"],
            verbose=args["verbose"],
            output_dir=args["stratification-output-dir"]
        )
    else
        logger = ConsoleLogger(stdout,Logging.Info)
        global_logger(logger)

        @info "MILK"
        flush(stdout)

        if isnothing(args["label"])
            label = replace(basename(args["input-path"]),r"\.csv(\.gz)?$" => "")
        else
            label = args["label"]
        end

        absolute_input_path = isabspath(args["input-path"]) ? args["input-path"] : joinpath(pwd(),args["input-path"])
        conversion_start = time()
        milk_input_path = prepare_milk_input(absolute_input_path,label)
        conversion_time = round(time()-conversion_start,digits=2)
        if args["convert-input-only"]
            @info "Conversion complete (--convert-input-only): $(milk_input_path)"
            return
        end
    
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

        recursion_start = time()
        i = 0 # Initial recursion

        full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"

        input_path = joinpath(invariant_args["output-dir"],"$(full_label).input.ids")
        symlink(milk_input_path,input_path)
        symlink(binary_path(milk_input_path),binary_path(input_path))
        symlink(group_meta_path(milk_input_path),group_meta_path(input_path))

        n = countlines(input_path)
    
        @info "Iteration: $i ($n objects)"
        flush(stdout)

        if n <= args["partition-size"]
            recursive_processing_direct_execution(
                representatives_path=input_path,
                iteration=i,
                label=label,
                invariant_args=invariant_args
            )
        else
            representatives_path,groups_path = recursive_processing_framework(
                input_path=input_path,
                label=full_label,
                invariant_args=invariant_args
            )

            n = countlines(representatives_path)
            while n > args["sample-size"]
                if count_rows(representatives_path) <= args["partition-size"]
                    recursive_processing_direct_execution(
                        representatives_path=representatives_path,
                        iteration=i+1,
                        label=label,
                        invariant_args=invariant_args
                    )
                    break
                end

                i += 1
                full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
                @info "Iteration: $i ($n objects)"
                flush(stdout)

                input_path = joinpath(invariant_args["output-dir"],"$(full_label).input.ids")
                symlink(representatives_path,input_path)
                symlink(binary_path(representatives_path),binary_path(input_path))
                symlink(group_meta_path(representatives_path),group_meta_path(input_path))

                if isfile(cache_meta_path(representatives_path))
                    symlink(cache_meta_path(representatives_path),cache_meta_path(input_path))
                    symlink(cache_binary_path(representatives_path),cache_binary_path(input_path))
                end

                representatives_path,groups_path = recursive_processing_framework(
                    input_path=input_path,
                    label=full_label,
                    invariant_args=invariant_args
                )
                n = countlines(representatives_path)
                @info "\t$n objects after recursive iteration."
                flush(stdout)
            end
        end

        @info "\nCompleted recursive downsampling procedure."
        flush(stdout)

        final_cleanup(invariant_args["output-dir"])
        recursion_time = round(time()-recursion_start,digits=2)

        reconstruction_start = time()
        if args["skip-reconstruction"]
            @info "Skipping cell hierarchy reconstruction."
        else
            @info "Reconstructing hierarchical graph..."
            hierarchical_reconstruction(args["output-dir"])
        end
        reconstruction_time = round(time()-reconstruction_start,digits=2)
        @info "Runtime (seconds): conversion $conversion_time; recursion $recursion_time; reconstruction $reconstruction_time; total $(round(conversion_time+recursion_time+reconstruction_time,digits=2))"
    
        @info "\tDone!"
        flush(stdout)
    end
end
