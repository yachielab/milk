
function instantitate_array_job_script(;job_scheduler,job_account,job_name,job_time,job_memory,num_cpus,stdout_dir,stderr_dir,environment_path)
    if job_scheduler == "sge"
        script_list = [
            "#!/bin/bash",
            "#\$ -S /bin/bash",
            "#\$ -N $(job_name)",
            "#\$ -pe def_slot $(num_cpus)",
            "#\$ -l s_vmem=$(job_memory)",
            "#\$ -o $(job_stdout_dir)",
            "#\$ -e $(job_stderr_dir)\n",
            "source $(environment_path)"
        ]
    elseif job_scheduler == "slurm"
        script_list = [
            "#!/bin/bash",
            "#SBATCH --account=$(job_account)",
            "#SBATCH --job-name$(job_name)",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=$(job_memory)",
            "#SBATCH --cpus-per-task=$(num_cpus)",
            "#SBATCH --time=$(job_time)",
            "#SBATCH --output=$(job_stdout_dir)",
            "#SBATCH --error=$(job_stderr_dir)\n",
            "source $(environment_path)"
        ]
    else
        throw
    end
    return script_list
end

function submit_stratification_batch_array_jobs()
    # pass
end
