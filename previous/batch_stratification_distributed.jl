# using ArgParse
# using StatsBase
# using Distributed
# using Logging
# using Glob

# @everywhere include("utils/utils.jl")
# @everywhere include("utils/file_handling.jl")
# @everywhere include("utils/pairwise_comparisons.jl")

function parse_arguments()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--input-dir","-i"
            arg_type = String
            required = true
            help = ""
        "--threshold","-t"
            arg_type = Float64
            default = nothing
            help = ""
        "--percentile","-p"
            arg_type = Float64
            default = nothing
            help = ""
        "--metric","-m"
            arg_type = String
            required = true
            help = ""
        "--cache-path","-c"
            arg_type = String
            default = nothing
            help = ""
        "--previous-groups-path","-G"
            arg_type = String
            default = nothing
            help = ""
        "--verbose","-v"
            action = :store_true 
        "--output-dir","-o"
            arg_type = String
            required = true
            help = ""
    end
    return parse_args(args)
end

function main()
    args = parse_arguments()

    @info "Starting distributed stratification process"
    @info "  $(nworkers())) workers"
    @info "  Input batch directory: $(args["input-dir"])"
    @info "  Threshold: $(args["threshold"])"
    @info "  Percentile: $(args["percentile"])"
    @info "  Metric: $(args["metric"])"
    @info "  Cache path: $(args["cache-path"])"
    @info "  Previous groups path: $(args["previous-groups-path"])"
    @info "  Output directory: $(args["output-dir"])"
    flush(stdout)

    n_comparisons = 0
    distance_function = map_distance_function(args["metric"])

    cache_dict = nothing
    if !isnothing(args["cache-path"]) && isfile(args["cache-path"])
        cache_dict = load_input_array_as_dictionary(args["cache-path"])
    end

    previous_groups = nothing
    if !isnothing(args["previous-groups-path"]) && isfile(args["previous-groups-path"])
        previous_groups_dict = load_groups_as_dictionary(args["previous-groups-path"])
        previous_groups = previous_groups_dict["groups"]
        previous_groups_dict = nothing
        GC.gc()
    end

    @info "Initiating the pmap distributed call"
    flush(stdout)
    pathlist = glob("*.csv",args["input-dir"])
    info_list = pmap(path -> stratification_process_distributed_execution(path,args["threshold"],args["percentile"],distance_function,cache_dict,previous_groups,args["output-dir"]), pathlist)

    if args["verbose"]
        for info in info_list
            @info info
        end
    end

end

main()
