# using ArgParse
# using StatsBase
# using Distributed
# using Logging
# using Glob

# @everywhere include("utils/utils.jl")
# @everywhere include("utils/file_handling.jl")
# @everywhere include("utils/pairwise_comparisons.jl")

function group_stratification(input_dir,threshold,percentile,metric,cache_path,
                              previous_groups_path,verbose,output_dir)

    @info "Starting distributed stratification process"
    @info "  $(nworkers()) workers"
    @info "  Input batch directory: $(input_dir)"
    @info "  Threshold: $(threshold)"
    @info "  Percentile: $(percentile)"
    @info "  Metric: $(metric)"
    @info "  Cache path: $(cache_path)"
    @info "  Previous groups path: $(previous_groups_path)"
    @info "  Output directory: $(output_dir)"
    flush(stdout)

    n_comparisons = 0
    distance_function = map_distance_function(metric)

    cache_dict = nothing
    if !isnothing(cache_path) && isfile(cache_path)
        cache_dict = load_input_array_as_dictionary(cache_path)
    end

    previous_groups = nothing
    if !isnothing(previous_groups_path) && isfile(previous_groups_path)
        previous_groups_dict = load_groups_as_dictionary(previous_groups_path)
        previous_groups = previous_groups_dict["groups"]
        previous_groups_dict = nothing
        GC.gc()
    end

    @info "Initiating the pmap distributed call"
    flush(stdout)
    pathlist = glob("*.csv",input_dir)
    info_list = pmap(path -> distributed_stratification_process(path,threshold,percentile,distance_function,cache_dict,previous_groups,output_dir), pathlist)

    if verbose
        for info in info_list
            @info info
        end
    end

end
