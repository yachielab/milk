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

    function julia_main()::Cint
        try
            args = parse_arguments()
            n_cpus = args["threads"]

            if nworkers() < n_cpus
                addprocs(n_cpus - nworkers())
            end

            for p in workers()
                Distributed.remotecall_eval(@__MODULE__, p, :(using Milk))
            end

            main()
        catch e
            Base.invokelatest(Base.display_error, e, catch_backtrace())
            return 1
        end
        return 0
    end
end
