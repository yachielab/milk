using Glob
include("file_handling.jl")

function instantitate_array_job_script(label,b,invariant_args)
    jobs_dir = joinpath(invariant_args["output-dir"],"jobs")
    stdout_dir = joinpath(jobs_dir,"stdout")
    stderr_dir = joinpath(jobs_dir,"stderr")
    mkpath(stdout_dir)
    mkpath(stderr_dir)
    if invariant_args["job-scheduler"] == "sge"
        job_script = """#!/bin/bash
        #\$ -S /bin/bash
        #\$ -N $(invariant_args["job-name"])
        #\$ -pe def_slot $(invariant_args["threads"])
        #\$ -l s_vmem=$(invariant_args["job-memory"])G
        #\$ -t 1:$(b)
        #\$ -o $(stdout_dir)
        #\$ -e $(stderr_dir)\n
        source $(invariant_args["environment-path"])\n
        ARRAY_INPUTS_PATH=\$1
        PERCENTILE=\$2
        METRIC=\$3
        THREADS=\$4
        OUTPUT_DIR=\$5\n
        BATCH_DIR=\$( sed -n \${SGE_TASK_ID}p \$ARRAY_INPUTS_PATH )\n
        """
    elseif invariant_args["job-scheduler"] == "slurm"
        job_script = """#!/bin/bash
        #SBATCH --account=$(invariant_args["job-account"])
        #SBATCH --job-name$(invariant_args["job-name"])
        #SBATCH --ntasks=1
        #SBATCH --nodes=1
        #SBATCH --mem=$(invariant_args["job-memory"])GB
        #SBATCH --cpus-per-task=$(invariant_args["threads"])
        #SBATCH --time=$(invariant_args["job-time"])
        #SBATCH --array=1-$(b)
        #SBATCH --output=$(stdout_dir)
        #SBATCH --error=$(stderr_dir)\n
        source $(invariant_args["environment-path"])\n
        ARRAY_INPUTS_PATH=\$1
        PERCENTILE=\$2
        METRIC=\$3
        THREADS=\$4
        OUTPUT_DIR=\$5\n
        BATCH_DIR=\$( sed -n \${SLURM_ARRAY_TASK_ID}p \$ARRAY_INPUTS_PATH )\
        """
    else
        throw
    end

    job_script *= """
    for PARTITIONED_INPUT_PATH in \$BATCH_DIR/*.csv.gz
    do
        LABEL=\$(basename -s ".csv.gz" \$PARTITIONED_INPUT_PATH)
        julia --threads \$THREADS $(dirname(@__DIR__))/threshold_based_group_stratification.jl \\
        -i \$PARTITIONED_INPUT_PATH \\
        -t \$PERCENTILE \\
        -l \$LABEL \\
        -m \$METRIC \\
        -o \$OUTPUT_DIR
    done
    """

    script_path = joinpath(jobs_dir,"$(label).array_jobs.sh")
    write_value_as_txt(job_script,script_path)
    return script_path
end

function submit_stratification_batch_array_jobs(;batches,files,label,partition_dir,invariant_args)
    # pass
    n = length(files)
    b = length(batches)
    array_inputs_path = joinpath(partition_dir,"$(label).batch_pathlist.txt")
    write_values_as_txt(batches,array_inputs_path)
    script_path = instantitate_array_job_script(label,b,invariant_args)
    command = [
        script_path,
        array_inputs_path,
        string(invariant_args["percentile"]),
        invariant_args["metric"],
        string(invariant_args["threads"]),
        partition_dir
    ]
    if invariant_args["job-scheduler"] == "sge"
        command = Cmd(vcat("qsub",command))
    elseif invariant_args["job-scheduler"] == "slurm"
        command = Cmd(vcat("sbatch",command))
    else
        throw
    end
    run(command)

    pattern = "$(label).*.groups.jsonl.gz"
    poll_for_job_completion(pattern,partition_dir,n)
    return
end

function poll_for_job_completion(pattern,partition_dir,n,polling_interval=30)
    @info "\tPolling for $n output files..."
    start_time = time()

    while length(glob(pattern,partition_dir)) < n
        sleep(polling_interval)
    end
    
    elapsed_time = round((time()-start_time)/60,digits=2)
    @info "\t\tFinished polling! (elapsed time: $elapsed_time min)"
    return
    end