function main()
    logger = ConsoleLogger(stdout,Logging.Info)
    global_logger(logger)

    args = parse_arguments()
    validate_args(args)

    if args["group-stratification-mode"]
        batch_group_stratification_hpc_mode(
            input_dir=args["stratification-input-dir"],
            parameters_path=args["stratification-parameters"]
        )
    else
        @info "MILK"
        flush(stdout)

        if isnothing(args["label"])
            label = replace(basename(args["input-path"]),r"\.csv(\.gz)?$" => "")
        else
            label = args["label"]
        end

        absolute_input_path = isabspath(args["input-path"]) ? args["input-path"] : joinpath(pwd(),args["input-path"])
        conversion_start = time()
        if isnothing(args["incremental-path"])
            milk_input_path = prepare_milk_input(absolute_input_path,label)
        else
            milk_input_path = prepare_incremental_milk_input(absolute_input_path,args["incremental-path"],label)
        end
        conversion_time = round(time()-conversion_start,digits=2)
        if args["convert-input-only"]
            @info "Conversion complete (--convert-input-only): $(milk_input_path)"
            return
        end

        #= Placement
        The input (query) objects are placed onto the reference tree in --reference-dir; the updated tree goes to --output-dir
        =#
        if !isnothing(args["reference-dir"])
            @info "Using $(nworkers()) worker(s) for distributed processing"
            log_args(args)
            reference_dir = abspath(args["reference-dir"])
            output_dir = abspath(args["output-dir"])
            query_label = replace(basename(milk_input_path),r"\.ids$" => "")

            reference_leaves_pathlist = sort(glob("$(label).leaves.*.ids",reference_dir))
            if isempty(reference_leaves_pathlist)
                error("No leaves ($(label).leaves.*.ids) in $(reference_dir); is it a MILK output directory for label '$(label)'?")
            end
            if isfile(leaves_path(reference_dir,label,query_label))
                error("'$(query_label)' is already placed in $(reference_dir); rename the query file.")
            end
            reference_leaf_ids = Set(reduce(vcat,[readlines(path) for path in reference_leaves_pathlist]))
            overlap = count(in(reference_leaf_ids),readlines(milk_input_path))
            if overlap > 0
                error("$(overlap) query ids are already leaves of the reference")
            end
            if isdir(output_dir)
                if args["force-overwrite"]
                    @warn "Output directory already exists! Overwriting."
                    rm(output_dir,recursive=true)
                else
                    error("Output directory ($(output_dir)) already exists! Exiting.")
                end
            end
            mkpath(output_dir)

            placement_start = time()
            placement(reference_dir,output_dir,label,milk_input_path,reference_leaves_pathlist;query_label=query_label,seed=args["seed"],partition_size=args["partition-size"])
            # leaves of the updated tree: the reference's leaves and the query (so it can be the next reference)
            for reference_leaves_path in reference_leaves_pathlist
                symlink_milk_files(reference_leaves_path,joinpath(output_dir,basename(reference_leaves_path)))
            end
            symlink_milk_files(milk_input_path,leaves_path(output_dir,label,query_label))
            placement_time = round(time()-placement_start,digits=2)

            reconstruction_start = time()
            if args["skip-reconstruction"]
                @info "Skipping cell hierarchy reconstruction."
            else
                @info "Reconstructing hierarchical graph..."
                hierarchical_reconstruction(output_dir)
            end
            reconstruction_time = round(time()-reconstruction_start,digits=2)
            @info "Runtime (seconds): conversion $conversion_time; placement $placement_time; reconstruction $reconstruction_time; total $(round(conversion_time+placement_time+reconstruction_time,digits=2))"
            @info "\tDone!"
            flush(stdout)
            return
        end
    
        #= Incremental MILK
        Initial tree is built from the earliest stage, later stages are placed onto it
        =#
        stages_path = nothing
        initial_stage = nothing
        if !isnothing(args["incremental-path"])
            stages_path = milk_input_path
            fields = split(readline(stages_path),'\t')
            initial_stage = parse(Int,fields[1])
            milk_input_path = String(fields[2])
            @info "Incremental MILK: initial tree from stage $(initial_stage) ($(milk_input_path))"
        end
    
        @info "Using $(nworkers()) worker(s) for distributed processing"
        log_args(args)
        output_dir = abspath(args["output-dir"])
        milk_output_dir = output_dir
        if !isnothing(stages_path)
            milk_output_dir = joinpath(output_dir,"stage_$(lpad(initial_stage,8,'0'))")
        end
        invariant_args = instantiate_invariant_args(args)
        invariant_args["output-dir"] = milk_output_dir # set once, before any use

        if isdir(args["output-dir"])
            if args["force-overwrite"]
                @warn "Output directory already exists! Overwriting."
                rm(args["output-dir"],recursive=true)
            else
                error("Output directory ($(args["output-dir"])) already exists! Exiting.")
            end
        end
        mkpath(milk_output_dir)

        recursion_start = time()
        i = 0 # Initial recursion

        full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"

        input_path = joinpath(milk_output_dir,"$(full_label).input.ids")
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

                input_path = joinpath(milk_output_dir,"$(full_label).input.ids")
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

        final_cleanup(milk_output_dir;keep_representatives=args["keep-intermediate-files"])
        recursion_time = round(time()-recursion_start,digits=2)

        # leaves: the input objects of iteration 0 (placement onto this tree reads them)
        symlink_milk_files(milk_input_path,leaves_path(output_dir,label,replace(basename(milk_input_path),r"\.ids$" => "")))
        if !isnothing(stages_path)
            error("Incremental MILK is being rewritten on top of placement; stages after $(initial_stage) are not placed yet.")
        end

        reconstruction_start = time()
        if args["skip-reconstruction"]
            @info "Skipping cell hierarchy reconstruction."
        else
            @info "Reconstructing hierarchical graph..."
            hierarchical_reconstruction(output_dir)
        end
        reconstruction_time = round(time()-reconstruction_start,digits=2)
        @info "Runtime (seconds): conversion $conversion_time; recursion $recursion_time; reconstruction $reconstruction_time; total $(round(conversion_time+recursion_time+reconstruction_time,digits=2))"
    
        @info "\tDone!"
        flush(stdout)
    end
end
