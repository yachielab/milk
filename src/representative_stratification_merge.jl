using ArgParse
using StatsBase
using Base.Threads
using Logging

include("utils/utils.jl")
include("utils/file_handling.jl")
include("utils/pairwise_comparisons.jl")
include("utils/merge_helpers.jl")

function parse_arguments()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--input-path","-i"
            arg_type = String
            required = true
            help = ""
        "--groups-path","-g"
            arg_type = String
            default = nothing
            help = ""
        "--threshold","-t"
            required = true
            arg_type = Float64
        "--previous-groups-path","-G"
            arg_type = String
            default = nothing
            help = ""
        "--cache-path","-c"
            arg_type = String
            default = nothing
            help = ""
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

    label = args["label"]
    distance_function = map_distance_function(args["metric"])

    concat_representatives_dict = load_input_array_as_dictionary(args["input-path"])
    concat_groups_dict = load_groups_as_dictionary(args["groups-path"])

    concat_groups = concat_groups_dict["groups"]
    n_comparisons = 0
    # threshold = minimum(concat_groups_dict["thresholds"])
    threshold = args["threshold"]

    n = sum(length(group) for group in values(concat_groups))

    cache_dict = nothing
    cache_label = "no"
    if !isnothing(args["cache-path"]) && isfile(args["cache-path"])
        cache_dict = load_input_array_as_dictionary(args["cache-path"])
        cache_label = "yes"
    end

    merged_representatives,merged_groups,n_comps = stratification(concat_representatives_dict,threshold,distance_function)
    optimized_dict,optimization_set,x = optimize_representatives(concat_representatives_dict,merged_groups,cache_dict,distance_function)
    n_comparisons += x

    optimized_dict,optimized_groups,distances_dict,specificity_dict,x = stratification_predefined_medoids(
        data_dict=concat_representatives_dict,
        representative_dict=optimized_dict,
        cache_dict=cache_dict,
        threshold=threshold,
        distance_function=distance_function
    )
    n_comparisons += x

    direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )

    optimized_groups = compile_previous_groupings(optimized_groups,concat_groups)

    G = length(optimized_dict)
    g = length(optimization_set)
    t = round(threshold,digits=4)

    N = sum(length(group) for group in values(optimized_groups))

    end_time = time()

    if args["verbose"]
        @info "\t[MERGE: $label] $G groups ($g groups optimized); $n objects; threshold: $t; cache: $cache_label; runtime: $(round(end_time-start_time,digits=2)) seconds"
        flush(stdout)
    end

    representatives_path = joinpath(args["output-dir"],"$label.representatives.csv")
    groups_path = joinpath(args["output-dir"],"$label.groups.jsonl.gz")
    write_dictionary_as_csv(optimized_dict,representatives_path)
    write_group_results(
        path=groups_path,
        label=label,
        stage="representative_stratification_merge",
        cache_label=cache_label,
        compiled_label="yes",
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
    rm(args["input-path"])
    rm(args["groups-path"])
end

main()
