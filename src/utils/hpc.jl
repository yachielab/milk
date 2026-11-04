module HPC

    using Glob
    using JSON
    using ..FileHandling: write_values_as_txt

    export submit_stratification_batch_array_jobs,
           read_stratification_parameters

    function instantiate_array_job_script(;label,n_batches,array_inputs_path,parameters_path,invariant_args)
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
            PARAMETERS_PATH=$(parameters_path)\n
            BATCH_DIR=\$( sed -n \${SGE_TASK_ID}p \$ARRAY_INPUTS_PATH )\n
            export JULIA_NUM_THREADS=$(invariant_args["threads"])\n
            """
        elseif invariant_args["job-scheduler"] == "slurm"
            job_script = """#!/bin/bash
            #SBATCH --account=$(invariant_args["job-account"])
            #SBATCH --job-name=$(invariant_args["job-name"])
            #SBATCH --ntasks=1
            #SBATCH --nodes=1
            #SBATCH --mem=$(invariant_args["job-memory"]*invariant_args["threads"])GB
            #SBATCH --cpus-per-task=$(invariant_args["threads"])
            #SBATCH --time=$(invariant_args["job-time"])
            #SBATCH --array=1-$(n_batches)
            #SBATCH --output=$(stdout_dir)/%A_%a.out
            #SBATCH --error=$(stderr_dir)/%A_%a.err\n
            source $(invariant_args["environment-path"])\n
            ARRAY_INPUTS_PATH=$(array_inputs_path)
            PARAMETERS_PATH=$(parameters_path)\n
            BATCH_DIR=\$( sed -n \${SLURM_ARRAY_TASK_ID}p \$ARRAY_INPUTS_PATH )\n
            export JULIA_NUM_THREADS=$(invariant_args["threads"])\n
            """
        else
            error("Unknown job scheduler: $(invariant_args["job-scheduler"])")
        end

        job_script *= """
        milk --group-stratification-mode \\
        --input-path \$BATCH_DIR \\
        --stratification-input-dir \$BATCH_DIR \\
        --stratification-parameters \$PARAMETERS_PATH \\
        -T $(invariant_args["threads"])
        """

        script_path = joinpath(jobs_dir,"$(label).array_jobs.sh")
        write_values_as_txt([job_script],script_path)
        return script_path
    end

    function submit_stratification_batch_array_jobs(;batches,files,threshold,label,partition_dir,invariant_args)
        n = length(files) # number of .flag files to wait for
        array_inputs_path = joinpath(partition_dir,"batch_dirlist.txt")
        write_values_as_txt(batches,array_inputs_path)
        parameters_path = joinpath(partition_dir,"stratification_parameters.json")
        write_stratification_parameters(parameters_path,threshold,partition_dir,invariant_args)
        script_path = instantiate_array_job_script(
            label=label,
            n_batches=length(batches),
            array_inputs_path=array_inputs_path,
            parameters_path=parameters_path,
            invariant_args=invariant_args
        )
        if invariant_args["job-scheduler"] == "sge"
            command = Cmd(["qsub",script_path])
        elseif invariant_args["job-scheduler"] == "slurm"
            # command = Cmd(vcat("sbatch",command))
            command = Cmd(["sbatch",script_path])
        else
            error("Unknown job scheduler: $(invariant_args["job-scheduler"])")
        end
        # run(command)
        run(pipeline(command,stdout=devnull,stderr=devnull))

        pattern = "$(label).*.flag"
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

    #=
    settings every array job needs for one iteration's group stratification, written once by the main
    process into the partition directory (instead of passing each one through the job script)
    =#
    function write_stratification_parameters(parameters_path,threshold,output_dir,invariant_args)
        parameters = Dict(
            "threshold"     => threshold,
            "percentile"    => invariant_args["percentile"],
            "metric"        => invariant_args["metric"],
            "linkage"       => invariant_args["linkage"],
            "cache-size"    => invariant_args["cache-size"],
            "size-weighted" => invariant_args["size-weighted"],
            "verbose"       => invariant_args["verbose"],
            "output-dir"    => output_dir
        )
        open(parameters_path,"w") do io
            JSON.print(io,parameters,4)
        end
    end

    function read_stratification_parameters(parameters_path)
        return JSON.parsefile(parameters_path)
    end
end