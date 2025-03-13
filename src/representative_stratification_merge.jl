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
        "--output-dir","-o"
            arg_type = String
            required = true
            help = ""
    end
    return parse_args(args)
end

function main()
    args = parse_arguments()

    label = args["label"]
    distance_function = map_distance_function(args["metric"])

    concat_representatives_dict = load_input_array_as_dictionary(args["input-path"])
    concat_groups_dict = load_groups_as_dictionary(args["groups-path"])
    concat_groups = concat_groups_dict["groups"]
    n_comparisons = sum(concat_groups_dict["comparisons"])
    threshold = maximum(concat_groups_dict["thresholds"])

    cache_dict = nothing
    cache_label = "no"
    if !isnothing(args["cache-path"]) && isfile(args["cache-path"])
        cache_dict = load_input_array_as_dictionary(args["cache-path"])
        cache_label = "yes"
    end

    merged_representatives,merged_groups,n_comps = stratification(concat_representatives_dict,threshold,distance_function)
    optimization_set = Set(id for (id,group) in merged_groups if length(group) > 1)
    optimized_dict,x = optimize_representatives(concat_representatives_dict,merged_groups,optimization_set,cache_dict,distance_function)
    n_comparisons += x

    optimized_groups,unmapped_dict,distances_dict,specificity_dict,x = stratification_predefined_medoids(
        data_dict=concat_representatives_dict,
        representative_dict=optimized_dict,
        threshold=threshold,
        distance_function=distance_function
    )
    n_comparisons += x

    if !isempty(intersect(keys(optimized_dict),keys(unmapped_dict)))
        error("Overlapping object IDs in (optimized) representative_dict and unmapped_dict!")
    end
    unmapped_set = Set(keys(unmapped_dict))
    merge!(optimized_dict,unmapped_dict)
    for id in keys(unmapped_dict)
        optimized_groups[id] = [id]
    end

    n = length(concat_representatives_dict)
    G = length(optimized_dict)
    g = length(optimization_set)
    u = length(unmapped_set)
    t = round(threshold,digits=4)
    N = sum(length(group) for group in values(optimized_groups))
    @info "\t[MERGE: $label] $n objects ($N total); threshold: $t; $G groups ($g groups optimized; $u unmapped); $n_comparisons comparisons (cache: $cache_label)"

    representatives_path = joinpath(args["output-dir"],"$label.representatives.csv.gz")
    groups_path = joinpath(args["output-dir"],"$label.jsonl.gz")
    write_dictionary_as_csv(optimized_dict,representatives_path)
    write_group_results(
        path=groups_path,
        n_groups=G,
        groups=optimized_groups,
        optimization_set=optimization_set,
        unmapped_set=unmapped_set,
        distances_dict=distances_dict,
        specificity_dict=specificity_dict,
        threshold=threshold,
        n_comparisons=n_comparisons
    )

end

main()
