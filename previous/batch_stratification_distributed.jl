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

# function distributed_stratification_process(path,threshold,perc,distance_function,cache_dict,previous_groups,output_dir)

#     label = replace(basename(path),".csv" => "")

#     start_time = time()

#     n_comparisons = 0
#     data_dict = load_input_array_as_dictionary(path)
#     D,index_map,threshold_,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function)
#     n_comparisons += x

#     if isnothing(threshold)
#         threshold_label = "computed"
#         threshold = threshold_
#     else
#         threshold_label = "predefined"
#     end

#     representative_dict,groups = stratification_precomputed(data_dict,D,index_map,threshold)
#     optimized_dict,optimization_set,x = optimize_representatives(data_dict,groups,cache_dict,distance_function)
#     n_comparisons += x
#     optimized_dict,optimized_groups,distances_dict,specificity_dict,x = stratification_predefined_medoids_precomputed(
#         data_dict=data_dict,
#         representative_dict=optimized_dict,
#         D=D,
#         index_map=index_map,
#         cache_dict=cache_dict,
#         threshold=threshold,
#         distance_function=distance_function
#     )
#     n_comparisons += x
#     direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )

#     if isnothing(previous_groups)
#         previous_groups_label = "no"
#     else
#         previous_groups_label = "yes"
#         optimized_groups = compile_previous_groupings(optimized_groups,previous_groups)
#     end

#     end_time = time()

#     G = length(optimized_dict)
#     g = length(optimization_set)
#     t = round(threshold,digits=4)
#     N = sum(length(group) for group in values(optimized_groups))
#     if isnothing(cache_dict)
#         cache_label = "no"
#     else
#         cache_label = "yes"
#     end

#     elapsed_time = round(end_time-start_time,digits=2)

#     info = "\t[$label] $G groups ($g groups optimized); $n objects ($N total); threshold: $t ($threshold_label); $n_comparisons comparisons ($threads thread(s); cache: $cache_label; previous_groups: $previous_groups_label); runtime: $(elapsed_time) seconds"

#     representatives_path = joinpath(output_dir,"$(label).representatives.csv")
#     groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
#     write_dictionary_as_csv(optimized_dict,representatives_path)
#     write_group_results(
#         path=groups_path,
#         label=label,
#         stage="batch_stratification_process",
#         cache_label=cache_label,
#         compiled_label=previous_groups_label,
#         n_input_objects=N,
#         n_groups=G,
#         groups=optimized_groups,
#         optimization_set=optimization_set,
#         direct_groupsize_dict=direct_groupsize_dict,
#         distances_dict=distances_dict,
#         specificity_dict=specificity_dict,
#         threshold=threshold,
#         n_comparisons=n_comparisons
#     )

#     return info
# end

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
    info_list = pmap(path -> distributed_stratification_process(path,args["threshold"],args["percentile"],distance_function,cache_dict,previous_groups,args["output-dir"]), pathlist)

    if args["verbose"]
        for info in info_list
            @info info
        end
    end

end

main()
