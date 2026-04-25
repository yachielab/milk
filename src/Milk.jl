module Milk

    using ArgParse
    using Logging
    using Distributed
    using CodecZlib
    using Distances
    using Glob
    using JSON
    using StatsBase
    using Base.Iterators
    using Random

    include("utils/file_handling.jl")
    include("utils/hpc.jl")
    include("utils/pairwise_comparisons.jl")
    include("utils/representative_optimization.jl")
    include("utils/group_stratification.jl")
    include("utils/group_aggregation.jl")
    include("utils/hierarchical_reconstruction.jl")
    include("utils/cli.jl")
    include("utils/utils.jl")

    using .FileHandling
    using .PairwiseComparisons
    using .RepresentativeOptimization
    using .GroupAggregation
    using .GroupStratification
    using .HierarchicalReconstruction
    using .CLI
    using .HPC
    using .Utils

    include("main.jl")
    export main

    function _ping_workers(n_cpus::Integer)
        workers_needed = max(n_cpus, 1) - nworkers()
        if workers_needed <= 0
            return 0
        end
    
        get!(ENV, "JULIA_WORKER_TIMEOUT", "60")
        for attempt in 1:3
            try
                addprocs(workers_needed)
                return 0
            catch e
                @warn "bulk addproc attempt $attempt failed"
            end
        end
        @warn "bulk addprocs failed; falling back to one-by-one worker startup"
    
        for i in 1:workers_needed
            success = false
            for attempt in 1:3
                try
                    newp = addprocs(1)
                    success = true
                    break
                catch e
                    @warn "addprocs(1) failed" worker_index=i attempt=attempt exception=(e, catch_backtrace())
                    sleep(2 * attempt)
                end
            end
            if !success
                error("Failed to start worker $i after retries")
            end
        end
    end

    function julia_main()::Cint
        try
            args = parse_arguments()
            n_cpus = args["threads"]

            _ping_workers(n_cpus)

            for p in workers()
                Distributed.remotecall_eval(Main, p, :(using Milk))
            end

            main()
        catch e
            Base.invokelatest(Base.display_error, e, catch_backtrace())
            return 1
        end
        return 0
    end
end
