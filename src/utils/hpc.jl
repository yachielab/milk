module HPC

    using Glob
    using ..FileHandling: write_values_as_txt

    export submit_stratification_batch_array_jobs

    function instantiate_array_job_script(;label,n_batches,array_inputs_path,threshold,cache_path,previous_groups_path,output_dir,invariant_args)
        jobs_dir = joinpath(invariant_args["output-dir"],"jobs")
        stdout_dir = joinpath(jobs_dir,"stdout")
        stderr_dir = joinpath(jobs_dir,"stderr")
        mkpath(stdout_dir)
        mkpath(stderr_dir)
        if invariant_args["job-scheduler"] == "sge"
            job_script = """#!/bin/bash
            #\$ -S /bin/bash
            #\$ -cwd
            #\$ -N $(invariant_args["job-name"])
            #\$ -pe def_slot $(invariant_args["threads"])
            #\$ -l s_vmem=$(invariant_args["job-memory"])G
            #\$ -t 1:$(n_batches)
            #\$ -o $(stdout_dir)
            #\$ -e $(stderr_dir)\n
            source $(invariant_args["environment-path"])\n
            ARRAY_INPUTS_PATH=$(array_inputs_path)
            THRESHOLD=$(threshold)
            PERCENTILE=$(invariant_args["percentile"])
            METRIC=$(invariant_args["metric"])
            CACHE_PATH=$(cache_path)
            PREVIOUS_GROUPS_PATH=$(previous_groups_path)
            OUTPUT_DIR=$(output_dir)\n
            BATCH_DIR=\$( sed -n \${SGE_TASK_ID}p \$ARRAY_INPUTS_PATH )\n
            export JULIA_NUM_THREADS=$(invariant_args["threads"])\n
            """
        elseif invariant_args["job-scheduler"] == "slurm"
            job_script = """#!/bin/bash
            #SBATCH --account=$(invariant_args["job-account"])
            #SBATCH --job-name=$(invariant_args["job-name"])
            #SBATCH --ntasks=1
            #SBATCH --nodes=1
            #SBATCH --mem=$(invariant_args["job-memory"])GB
            #SBATCH --cpus-per-task=$(invariant_args["threads"])
            #SBATCH --time=$(invariant_args["job-time"])
            #SBATCH --array=1-$(n_batches)
            #SBATCH --output=$(stdout_dir)/%A_%a.out
            #SBATCH --error=$(stderr_dir)/%A_%a.err\n
            source $(invariant_args["environment-path"])\n
            ARRAY_INPUTS_PATH=$(array_inputs_path)
            THRESHOLD=$(threshold)
            PERCENTILE=$(invariant_args["percentile"])
            METRIC=$(invariant_args["metric"])
            CACHE_PATH=$(cache_path)
            PREVIOUS_GROUPS_PATH=$(previous_groups_path)
            OUTPUT_DIR=$(output_dir)\n
            BATCH_DIR=\$( sed -n \${SLURM_ARRAY_TASK_ID}p \$ARRAY_INPUTS_PATH )\n
            export JULIA_NUM_THREADS=$(invariant_args["threads"])\n
            """
        else
            throw
        end

        job_script *= """
        milk --group-stratification-mode \\
        --input-path \$BATCH_DIR \\
        --stratification-input-dir \$BATCH_DIR \\
        --stratification-threshold \$THRESHOLD \\
        --stratification-percentile \$PERCENTILE \\
        --stratification-metric \$METRIC \\
        --stratification-cache-path \$CACHE_PATH \\
        --stratification-previous-groups-path \$PREVIOUS_GROUPS_PATH \\
        -T $(invariant_args["threads"]) \\
        --stratification-output-dir \$OUTPUT_DIR
        """

        script_path = joinpath(jobs_dir,"$(label).array_jobs.sh")
        write_values_as_txt([job_script],script_path)
        return script_path
    end

    function submit_stratification_batch_array_jobs(;batches,files,threshold,cache_path,previous_groups_path,label,partition_dir,invariant_args)
        n = length(files)
        array_inputs_path = joinpath(partition_dir,"batch_dirlist.txt")
        write_values_as_txt(batches,array_inputs_path)
        if isnothing(cache_path)
            cache_path = "nothing"
        end
        if isnothing(previous_groups_path)
            previous_groups_path = "nothing"
        end
        script_path = instantiate_array_job_script(
            label=label,
            n_batches=length(batches),
            array_inputs_path=array_inputs_path,
            threshold=threshold,
            cache_path=cache_path,
            previous_groups_path=previous_groups_path,
            output_dir=partition_dir,
            invariant_args=invariant_args
        )
        if invariant_args["job-scheduler"] == "sge"
            command = Cmd(["qsub",script_path])
        elseif invariant_args["job-scheduler"] == "slurm"
            # command = Cmd(vcat("sbatch",command))
            command = Cmd(["sbatch",script_path])
        else
            throw
        end
        # run(command)
        run(pipeline(command,stdout=devnull,stderr=devnull))

        pattern = "$(label).*.groups.jsonl.gz"
        poll_for_job_completion(pattern,partition_dir,n)
        return
    end

    function poll_for_job_completion(pattern,partition_dir,n,polling_interval=30)
        @info "\tPolling for $n output files..."
        flush(stdout)
        start_time = time()

        while length(glob(pattern,partition_dir)) < n
            sleep(polling_interval)
        end
        
        elapsed_time = round((time()-start_time)/60,digits=2)
        @info "\t\tFinished polling! (elapsed time: $elapsed_time min)"
        flush(stdout)
        return
    end
end