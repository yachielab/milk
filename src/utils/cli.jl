module CLI
    
    using Logging
    using ArgParse

    export instantiate_invariant_args,log_args,parse_arguments,validate_args

    function parse_arguments()
        args = ArgParseSettings(description="A command-line tool to capture hierarchical relationships at scale.")
        add_arg_group(args, "Main arguments")
        @add_arg_table args begin
            "--input-path","-i" 
                arg_type = String
                required = true
                help = "Uncompressed CSV file. First column corresponds to object IDs (e.g., cell barcode). No header"
            "--sample-size","-n"
                arg_type = Int
                default = 1
                help = "Sample size threshold for the number of groups required to stop the recursive grouping process. Recursion stops when the number of (representative) objects is less than or equal to this value."
            "--cache-size-limit","-c"
                arg_type = Int
                default = 50000
                help = "How many (representative) objects to cache for medoid optimization in subsequent recursions."
            "--partition-size","-p"
                arg_type = Int
                default = 10000
                help = "How many objects per partitioned file."
            "--merge-threshold","-M"
                arg_type = Int
                default = 100000
                help = "Threshold of when to carry out an additional merging stratification process (if exceeded, basic concatenation to join partitioned groups)"
            "--metric","-m"
                arg_type = String
                default = "euclidean"
                help = "Distance metric for pairwise comparisons. Permissible values: {'cosine','euclidean','manhattan','hamming','jaccard','correlation'})"
            "--label","-l"
                arg_type = String
                default = nothing
                help = "Specify string label for file names."
            "--percentile","-t"
                arg_type = Float64
                default = 1.0
                help = "Threshold for stratification process to group cells."
            "--threads","-T"
                arg_type = Int
                default = 1
                help = "Number of distributed processes for computing pairwise comparisons. (TODO change to n_cpus)"
            "--seed","-r"
                arg_type = Int
                default = 21
                help = "Random seed."
            "--verbose"
                action = :store_true
                help = "Amount of information written to standard out/err."
            "--force-overwrite"
                action = :store_true
                help = "Overwrite output directory if it already exists."
            "--skip-reconstruction"
                action = :store_true
                help = "Only perform hierarchical grouping/downsampling phase and NOT the additional step of reconstructing vertices and edges table for graph structure."
            "--output-dir","-o"
                arg_type = String
                default = "./milk.out"
                help = "Directory to write intermediate and output files."
        end

        add_arg_group(args,"[High-Performance Computing mode] Cluster computing for large-scale MILK runs")
        @add_arg_table args begin
            "--hpc-mode"
                action = :store_true
                help = "Activate HPC mode for MILK execution. If this flag is not specified, all [HPC mode] arguments will be ignored."
            "--batch-size","-b"
                arg_type = Int
                default = 50
                help = "[HPC mode] The number of partitioned input files to be assigned per job."
            "--job-scheduler"
                arg_type = String
                default = nothing
                help = "[HPC mode] Job scheduler to distribute partitioned batches of stratification processes. Supported schedulers: {'slurm','sge'}"
            "--job-account"
                arg_type = String
                default = nothing
                help = "[HPC mode] Account name associated with account (only for SLURM job scheduler)."
            "--job-name"
                arg_type = String
                default = "milk"
                help = "[HPC mode] Job names"
            "--job-memory"
                arg_type = Int
                default = 4
                help = "[HPC mode] Memory of jobs (in gigabytes)."
            "--job-time"
                arg_type = String
                default = "24:00:00"
                help = "[HPC mode] Allotted time for jobs (only for SLURM job scheduler)."
            "--environment-path"
                arg_type = String
                default = nothing
                help = "[HPC mode] Bash script to set up environment in submitted jobs."
        end

        add_arg_group(args,"[for INTERNAL calls only] distributed batch group stratification")
        @add_arg_table args begin
            "--group-stratification-mode"
                action = :store_true
            "--stratification-input-dir"
                arg_type = String
            "--stratification-threshold"
                arg_type = Float64
            "--stratification-percentile"
                arg_type = Float64
                default = nothing
            "--stratification-metric"
                arg_type = String
            "--stratification-cache-path"
                arg_type = String
                default = nothing
            "--stratification-previous-groups-path"
                arg_type = String
                default = nothing
            "--stratification-output-dir"
                arg_type = String
        end

        return parse_args(args)
    end

    function instantiate_invariant_args(args)
        invariant_args = Dict(
            "batch-size"       => args["batch-size"],
            "sample-size"      => args["sample-size"],
            "cache-size-limit" => args["cache-size-limit"],
            "partition-size"   => args["partition-size"],
            "merge-threshold"  => args["merge-threshold"],
            "percentile"       => args["percentile"],
            "metric"           => args["metric"],
            "seed"             => args["seed"],
            "threads"          => args["threads"],
            "job-scheduler"    => args["job-scheduler"],
            "job-account"      => args["job-account"],
            "job-name"         => args["job-name"],
            "job-time"         => args["job-time"],
            "job-memory"       => args["job-memory"],
            "environment-path" => isnothing(args["environment-path"]) ? nothing : abspath(args["environment-path"]),
            "verbose"          => args["verbose"],
            "output-dir"       => abspath(args["output-dir"]),
            "hpc-mode"         => args["hpc-mode"]
        )
        return invariant_args
    end

    function log_args(args)
        @info "Parsed Arguments:"
        for (key,value) in args
            if startswith(key,"stratification-") continue end
            @info "  $(key): $(value)"
        end
        @info "=========="
        flush(stdout)
    end 

    function validate_args(args)
        if args["hpc-mode"]
            supported_job_schedulers = Set{String}(["slurm","sge"])
            if !(args["job-scheduler"] in supported_job_schedulers)
                error("$(args["job_scheduler"]) is either invalid or not supported.")
            end
            if args["job-scheduler"] == "slurm"
                if isnothing(args["job-account"])
                    error("A job account name is required for SLURM execution.")
                end
            end
            if (isnothing(args["environment-path"]) || !isfile(args["environment-path"]))
                error("A bash script to initialize environmental variables in submitted jobs is required!")
            end
        end
        if args["group-stratification-mode"]
            if (isnothing(args["stratification-input-dir"]) || !isdir(args["stratification-input-dir"]))
                error("Invalid '--stratification-input-dir' specified ($(args["stratification-input-dir"])).")
            end
            if (isnothing(args["stratification-output-dir"]) || !isdir(args["stratification-output-dir"]))
                error("Invalid '--stratification-output-dir' specified ($(args["stratification-output-dir"])).")
            end
        end
    end


end
