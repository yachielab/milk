using ArgParse
using StatsBase
using Base.Threads
using Logging

include("utils/utils.jl")
include("utils/file_handling.jl")
include("utils/pairwise_comparisons.jl")

function parse_arguments()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--input-path","-i"
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
        "--threshold","-t"
            arg_type = Float64
            default = nothing
            help = ""
        "--percentile","-p"
            required = true
            arg_type = Float64
            help = "Ignored if `-t` threshold argument is specified."
        "--label","-l"
            arg_type = String
            default = nothing
            help = ""
        "--metric","-m"
            arg_type = String
            required = true
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

    start_time = time()

    n_comparisons = 0
    threads = Threads.nthreads()
    threshold = args["threshold"]
    distance_function = map_distance_function(args["metric"])

    if isnothing(args["label"])
        label = replace(basename(args["input-path"]),".csv.gz" => "")
    else
        label = args["label"]
    end

    data_dict = load_input_array_as_dictionary(args["input-path"])
    n = length(data_dict)

    # TODO: error when there is only one sequence in the partitioned input file
    D,index_map,threshold,x = pairwise_distances(data_dict,threads,args["percentile"],distance_function)
    flush(stdout)
    n_comparisons += x

    if !isnothing(args["threshold"])
        threshold_label = "predefined"
        threshold = args["threshold"]
    else
        threshold_label = "computed"
    end

    representative_dict,groups = stratification_precomputed(data_dict,D,index_map,threshold)

    cache_dict = nothing
    cache_label = "no"
    if !isnothing(args["cache-path"]) && isfile(args["cache-path"])
        cache_dict = load_input_array_as_dictionary(args["cache-path"])
        cache_label = "yes"
    end

    previous_groups_label = "no"
    if !isnothing(args["previous-groups-path"]) && isfile(args["previous-groups-path"])
        previous_groups_dict = load_groups_as_dictionary(args["previous-groups-path"])
        previous_groups = previous_groups_dict["groups"]
        groups = compile_previous_groupings(groups,previous_groups)
        previous_groups_label = "yes"
    end

    optimized_dict,optimization_set,x = optimize_representatives(data_dict,groups,cache_dict,distance_function)
    n_comparisons += x

    optimized_dict,optimized_groups,distances_dict,specificity_dict,x = stratification_predefined_medoids_precomputed(
        data_dict=data_dict,
        representative_dict=optimized_dict,
        D=D,
        index_map=index_map,
        cache_dict=cache_dict,
        threshold=threshold,
        distance_function=distance_function
    )

    direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )

    if previous_groups_label == "yes"
        optimized_groups = compile_previous_groupings(optimized_groups,previous_groups)
    end

    n_comparisons += x

    G = length(optimized_dict)
    g = length(optimization_set)
    t = round(threshold,digits=4)
    N = sum(length(group) for group in values(optimized_groups))

    end_time = time()

    if args["verbose"]
        @info "\t[$label] $G groups ($g groups optimized); $n objects ($N total); threshold: $t ($threshold_label); $n_comparisons comparisons ($threads thread(s); cache: $cache_label; previous_groups: $previous_groups_label); runtime: $(round(end_time-start_time,digits=2)) seconds"
        flush(stdout)
    end

    representatives_path = joinpath(args["output-dir"],"$(args["label"]).representatives.csv.gz")
    groups_path = joinpath(args["output-dir"],"$(args["label"]).groups.jsonl.gz")
    write_dictionary_as_csv(optimized_dict,representatives_path)
    write_group_results(
        path=groups_path,
        label=label,
        stage="threshold_based_group_stratification",
        cache_label=cache_label,
        compiled_label=previous_groups_label,
        n_input_objects=N,
        n_groups=G,
        groups=optimized_groups,
        optimization_set=optimization_set,
        direct_groupsize_dict=direct_groupsize_dict,
        distances_dict=distances_dict,
        specificity_dict=specificity_dict,
        threshold=threshold,
        n_comparisons=n_comparisons
    )
end

main()
