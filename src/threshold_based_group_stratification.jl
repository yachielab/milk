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
        "--previous-groups-path","-p"
            arg_type = String
            default = nothing
            help = ""
        "--cache-path","-c"
            arg_type = String
            default = nothing
            help = ""
        "--percentile","-t"
            required = true
            default = 1.0
            arg_type = Float64
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

    n_comparisons = 0
    threads = Threads.nthreads()
    distance_function = map_distance_function(args["metric"])

    if isnothing(args["label"])
        label = replace(basename(args["input-path"]),".csv.gz" => "")
    else
        label = args["label"]
    end

    ids,profiles = load_input_array_as_lists(args["input-path"])
    index_map = Dict( ids[i] => i for i in eachindex(ids) )
    n = length(ids)

    D,threshold,n_pairwise_comps = pairwise_distances(profiles,threads,args["percentile"],distance_function)
    n_comparisons += n_pairwise_comps
    representative_dict,groups = group_stratification_precomputed_distances(
        ids=ids,
        vecs=profiles,
        D=D,
        index_map=index_map,
        threshold=threshold
    )
    n_groups = length(groups)

    optimization_set = Set(id for (id,group) in groups if length(group) > 1)

    cache_dict = nothing
    cache_label = "no"
    if !isnothing(args["cache-path"]) && isfile(args["cache-path"])
        cache_dict = load_input_array_as_dictionary(args["cache-path"])
        cache_label = "yes"
    end

    optimized_dict,optimized_groups,unmapped_set,distances_dict,specificity_dict,n_opt_comps = optimization_process(
        ids=ids,
        vecs=profiles,
        D=D,
        index_map=index_map,
        groups=groups,
        optimization_set=optimization_set,
        cache_dict=cache_dict,
        threshold=threshold,
        distance_function=distance_function
    )
    
    if !isnothing(args["previous-groups-path"]) && isfile(args["previous-groups-path"])
        previous_groups = load_groups_as_dictionary(args["previous-groups-path"])
        compiled_groups,optimization_set = compile_previous_groupings(optimized_groups,previous_groups)
        optimized_groups = compiled_groups
        previous_groups_label = "yes"
    else
        previous_groups_label = "no"
    end

    n_groups_optimized = length(optimization_set)
    n_unmapped = length(unmapped_set)
    rounded_threshold = round(threshold,digits=4)
    n_comparisons += n_opt_comps
    N = sum(length(group) for group in values(optimized_groups))
    @info "\t[$label] $n objects ($N total); threshold: $rounded_threshold; $n_groups groups ($n_groups_optimized groups optimized; $n_unmapped unmapped); $n_comparisons comparisons ($threads thread(s); cache: $cache_label; previous groups: $previous_groups_label)"

    representatives_path = joinpath(args["output-dir"],"$(args["label"]).representatives.csv.gz")
    groups_path = joinpath(args["output-dir"],"$(args["label"]).groups.jsonl.gz")
    write_dictionary_as_csv(optimized_dict,representatives_path)
    write_group_results(
        path=groups_path,
        n_groups=n_groups,
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
