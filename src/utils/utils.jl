# using GZip
# using Glob
# using CodecZlib
# using Random
# using StatsBase
# using Distances
# using Base.Iterators
# using Base.Threads
# using Logging
# using Distributed

# include("hpc.jl")
# include("file_handling.jl")
# @everywhere include("pairwise_comparisons.jl")

@everywhere function map_object_precomputed(i,representative_dict,D,index_map,threshold)
    """
    Map each candidate (denoted by its index, i) to the representative it is closest to
    (This returns nothing if there are no representatives it is similar to, as defined by
    the threshold)
    """
    closest_id = nothing
    closest_dist = Inf
    specificity = 0
    for (representative_id,representative_vec) in representative_dict
        j = index_map[representative_id]
        @inbounds d = D[i,j]
        if d <= threshold
            specificity += 1
            if d < closest_dist
                closest_id = representative_id
                closest_dist = d
            end
        end
    end
    return closest_id,closest_dist,specificity
end

@everywhere function stratification_precomputed(data_dict,D,index_map,threshold)
    representative_dict = Dict{String,Vector{Float32}}()
    groups = Dict{String,Vector{String}}()
    for (candidate_id,candidate_vec) in data_dict
        i = index_map[candidate_id]
        closest_id,closest_dist,_ = map_object_precomputed(i,representative_dict,D,index_map,threshold)
        if isnothing(closest_id)
            representative_dict[candidate_id] = candidate_vec
            groups[candidate_id] = [candidate_id]
        else
            push!(groups[closest_id],candidate_id)
        end
    end
    return representative_dict,groups
end

@everywhere function map_object(vec,representative_dict,threshold,distance_function)
    closest_id = nothing
    closest_dist = Inf
    specificity = 0
    n_comparisons = 0
    for (representative_id,representative_vec) in representative_dict
        d = distance_function(vec,representative_vec)
        n_comparisons += 1
        if d <= threshold
            specificity += 1
            if d < closest_dist
                closest_id = representative_id
                closest_dist = d
            end
        end
    end
    return closest_id,closest_dist,n_comparisons,specificity
end

@everywhere function stratification_predefined_medoids_precomputed(;data_dict,representative_dict,D,index_map,cache_dict,threshold,distance_function)
    groups = Dict( id => [id] for id in keys(representative_dict) )
    distances_dict = Dict(id => Vector{Float16}() for id in keys(representative_dict))
    specificity_dict = Dict(id => Vector{Int32}() for id in keys(representative_dict))
    for (candidate_id,candidate_vec) in data_dict
        if haskey(representative_dict,candidate_id) continue end
        i = index_map[candidate_id]
        closest_id,closest_dist,specificity = map_object_precomputed(i,representative_dict,D,index_map,threshold)
        if isnothing(closest_id)
            """
            New groups can still be formed if candidates are not similar to any of the optimized medoids, where
            similarity is determined by the threshold.
            """
            representative_dict[candidate_id] = candidate_vec
            groups[candidate_id] = [candidate_id]
            distances_dict[candidate_id] = Vector{Float16}()
            specificity_dict[candidate_id] = Vector{Int32}()
        else
            push!(groups[closest_id],candidate_id)
            push!(distances_dict[closest_id],closest_dist)
            push!(specificity_dict[closest_id],specificity)
        end
    end
    x = 0
    if !isnothing(cache_dict)
        for (cache_id,cache_vec) in cache_dict
            if haskey(index_map,cache_id) continue end
            closest_id,closest_dist,n_comps,specificity = map_object(cache_vec,representative_dict,threshold,distance_function)
            x += n_comps
            if !isnothing(closest_id)
                push!(distances_dict[closest_id],closest_dist)
                push!(specificity_dict[closest_id],specificity)
            end
        end
    end
    return representative_dict,groups,distances_dict,specificity_dict,x
end

@everywhere function compute_centroid(medoid_dict)
    return vec(mean(reduce(hcat,values(medoid_dict)),dims=2))
end

@everywhere function optimize_representatives(data_dict,groups,cache_dict,distance_function)
    optimized_dict = Dict{String,Vector{Float32}}()
    optimization_set = Set{String}()
    x = 0
    for (representative_id,group) in groups
        if length(group) == 1
            optimized_dict[representative_id] = data_dict[representative_id]
        else
            medoid_dict,cache_set = get_medoid_dictionary(group,data_dict,cache_dict)
            # if length(medoid_dict) < 2
            #     optimized_dict[representative_id] = data_dict[representative_id]
            # else
            #     centroid_vec = compute_centroid(medoid_dict)
            #     medoid_id,medoid_vec,n_centroid_comps = determine_medoid(medoid_dict,centroid_vec,cache_set,distance_function)
            #     x += n_centroid_comps
            #     optimized_dict[medoid_id] = medoid_vec
            #     push!(optimization_set,medoid_id)
            # end
            centroid_vec = compute_centroid(medoid_dict)
            medoid_id,medoid_vec,n_centroid_comps = determine_medoid(medoid_dict,centroid_vec,cache_set,distance_function)
            x += n_centroid_comps
            optimized_dict[medoid_id] = medoid_vec
            push!(optimization_set,medoid_id)
        end
    end
    return optimized_dict,optimization_set,x
end

function get_medoid_dictionary(group,representative_dict,cache_dict)
    medoid_dict = Dict{String,Vector{Float32}}() # this is for a single group
    cache_set = Set{String}()
    for sample_id in group
        if haskey(representative_dict,sample_id)
            medoid_dict[sample_id] = representative_dict[sample_id]
        elseif !isnothing(cache_dict) && haskey(cache_dict,sample_id)
            medoid_dict[sample_id] = cache_dict[sample_id]
            push!(cache_set,sample_id)
        end
    end
    return medoid_dict,cache_set
end

@everywhere function determine_medoid(medoid_dict,centroid_vec,cache_set,distance_function)
    medoid_id = nothing
    medoid_vec = nothing
    medoid_dist = Inf
    n_comparisons = 0
    for (sample_id,vec) in medoid_dict
        if sample_id in cache_set continue end
        dist = distance_function(centroid_vec,vec)
        n_comparisons += 1
        if dist < medoid_dist
            medoid_id   = sample_id
            medoid_vec  = vec
            medoid_dist = dist
        end
    end
    return medoid_id,medoid_vec,n_comparisons
end

function merge_groups(merged_groups,concat_groups)
    groups = Dict{String,Vector{String}}()
    for (representative_id,group) in merged_groups
        groups[representative_id] = Vector{String}()
        for id in group
            append!(groups[representative_id],concat_groups[id])
        end
    end
    return groups
end

function compile_previous_groupings(groups,previous_groups)
    compiled_groups = Dict{String,Vector{String}}()
    for (representative_id,group) in groups
        """
        Case when the previous and current groupings are exact same.
        No optimization of the medoid is required.
        """
        if (haskey(previous_groups,representative_id)) && (length(group) == 1)
            if group[1] != representative_id
                error("Unexpected state: group[1] ($(group[1])) is not the representative_id ($representative_id)")
            end
            compiled_groups[representative_id] = previous_groups[representative_id]
        else
            compiled_groups[representative_id] = Vector{String}()
            for id in group
                if !haskey(previous_groups, id)
                    error("Missing previous group for id: $id")
                end
                append!(compiled_groups[representative_id],previous_groups[id])
            end
        end
    end
    return compiled_groups
end


function attempt_to_merge_partitioned_results(;concat_representatives_path,concat_groups_path,p,m,threshold,label,cache_path,invariant_args)
    representatives_path = joinpath(invariant_args["output-dir"],"$(label).representatives.csv")
    groups_path = joinpath(invariant_args["output-dir"],"$(label).groups.jsonl.gz")
    if (p > 1) && (m <= invariant_args["merge-threshold"])
        # @info "\tMerging: number of concatenated groups, $m <= merge threshold ($(invariant_args["merge-threshold"]))"
        # command = [
        #     "julia",
        #     joinpath(dirname(@__DIR__),"representative_stratification_merge.jl"),
        #     "-i", concat_representatives_path,
        #     "-g", concat_groups_path,
        #     "-t", string(threshold),
        #     "-l", label,
        #     "-m", invariant_args["metric"],
        #     "--verbose",
        #     "-o", invariant_args["output-dir"]
        # ]
        group_aggregation_by_stratification(
            input_path=concat_representatives_path,
            groups_path=concat_groups_path,
            threshold=threshold,
            cache_path=cache_path,
            label=label,
            metric=invariant_args["metric"],
            output_dir=invariant_args["output-dir"]
        )
        command = Cmd(command)
        run(command)
    else
        mv(concat_representatives_path,representatives_path)
        mv(concat_groups_path,groups_path)
    end
    return representatives_path,groups_path
end


@everywhere function distributed_stratification_process(path,threshold,perc,distance_function,cache_dict,previous_groups,output_dir)

    label = replace(basename(path),".csv" => "")

    @info "$label"
    flush(stdout)
    # @info "$(nworkers())) workers"

    start_time = time()

    n_comparisons = 0
    data_dict = load_input_array_as_dictionary(path)
    D,index_map,threshold_,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function)
    n_comparisons += x

    if isnothing(threshold)
        threshold_label = "computed"
        threshold = threshold_
    else
        threshold_label = "predefined"
    end

    representative_dict,groups = stratification_precomputed(data_dict,D,index_map,threshold)
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
    n_comparisons += x
    direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )

    if isnothing(previous_groups)
        previous_groups_label = "no"
    else
        previous_groups_label = "yes"
         optimized_groups = compile_previous_groupings(optimized_groups,previous_groups)
    end

    end_time = time()

    n = length(data_dict)
    G = length(optimized_dict)
    g = length(optimization_set)
    t = round(threshold,digits=4)
    N = sum(length(group) for group in values(optimized_groups))
    if isnothing(cache_dict)
        cache_label = "no"
    else
        cache_label = "yes"
    end

    elapsed_time = round(end_time-start_time,digits=2)

    info = "\t[$label] $G groups ($g groups optimized); $n objects ($N total); threshold: $t ($threshold_label); $n_comparisons comparisons ($(nworkers())) CPU(s); cache: $cache_label; previous_groups: $previous_groups_label); runtime: $(elapsed_time) seconds"

    representatives_path = joinpath(output_dir,"$(label).representatives.csv")
    groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
    write_dictionary_as_csv(optimized_dict,representatives_path)
    write_group_results(
        path=groups_path,
        label=label,
        stage="batch_stratification_process",
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
    rm(path)
    return info
end

function stratify_partitioned_populations_into_groups(;batches,files,threshold,cache_path,previous_groups_path,label,partition_dir,invariant_args)
    b = length(batches)
    perc = invariant_args["percentile"]
    distance_function = map_distance_function(invariant_args["metric"])
    cache_dict = nothing
    if !isnothing(cache_path) && isfile(cache_path)
        cache_dict = load_input_array_as_dictionary(cache_path)
    end
    previous_groups = nothing
    if !isnothing(previous_groups_path) && isfile(previous_groups_path)
        previous_groups_dict = load_groups_as_dictionary(previous_groups_path)
        previous_groups = previous_groups_dict["groups"]
    end
    if b == 1
        batch_dir = batches[1]
        pattern = "*.csv"
        pathlist = sort(glob(pattern,batch_dir))
        info_list = pmap(path -> distributed_stratification_process(path,threshold,perc,distance_function,cache_dict,previous_groups,partition_dir), pathlist)
        @everywhere GC.gc() 
        if invariant_args["verbose"]
            for info in info_list
                @info info
            end
        end
    else
        @info "\tSubmitting $b batches as an array job..."
        submit_stratification_batch_array_jobs(
            batches=batches,
            files=files,
            threshold=threshold,
            cache_path=cache_path,
            previous_groups_path=previous_groups_path,
            label=label,
            partition_dir=partition_dir,
            invariant_args=invariant_args
        )
    end
    return
end

function attempt_to_cache_file(path,n,invariant_args)
    if n <= invariant_args["cache-size-limit"]
        @info "\tCaching $n objects from the following path: $path" 
        flush(stdout)
        return path
    else
        return nothing
    end
end

function process(;input_path,label,cache_path,previous_groups_path,invariant_args)

    # @info "\tPartitioning input file..."
    partition_dir,files,batches = partition_and_batch_input_files(input_path,label,invariant_args)
    n_batches = length(batches)
    @info "\tInput file partitioned into $(length(files)) files ($(length(batches)) batches)"
    flush(stdout)

    @info "\tDetermining threshold: randomly sampling files..." 
    flush(stdout)
    if length(files) > 1
        if length(files) > invariant_args["batch-size"]
            Random.seed!(invariant_args["seed"])
            sampled_files = sample(files,invariant_args["batch-size"])
        else
            sampled_files = files
        end
        start_time = time()
        threshold = determine_percentile_threshold_distributed(sampled_files,invariant_args)
        elapsed_time = round(time()-start_time,digits=2)
        @info "\tGlobally-determined percentile threshold: $(threshold) ($(length(sampled_files)) files; $elapsed_time seconds)"
        flush(stdout)
    else
        threshold = nothing
    end

    # @info "\tstratifying partitioned input files..."
    stratify_partitioned_populations_into_groups(
        batches=batches,
        files=files,
        threshold=threshold,
        cache_path=cache_path,
        previous_groups_path=previous_groups_path,
        label=label,
        partition_dir=partition_dir,
        invariant_args=invariant_args
    )

    concat_representatives_path,concat_groups_path,m = concatenate_partitioned_results(label,partition_dir,invariant_args)
    @info "\t$m after concatenation across partitioned results."
    flush(stdout)

    representatives_path,groups_path = attempt_to_merge_partitioned_results(
        concat_representatives_path=concat_representatives_path,
        concat_groups_path=concat_groups_path,
        p=length(files),
        m=m,
        threshold=threshold,
        label=label,
        cache_path=cache_path,
        invariant_args=invariant_args
    )
    exclusion_set = Set{String}([representatives_path])
    if !isnothing(cache_path) && isfile(cache_path)
        push!(exclusion_set,cache_path)
    end
    clean_directory(invariant_args["output-dir"],partition_dir,exclusion_set)
    return representatives_path,groups_path
end


function direct_recursive_processing(;representatives_path,iteration,label,cache_path,previous_groups_path,invariant_args)
    # carrying out remaining recursions in a single function rather than the framework (overhead is high for single input file)
    @info "========== Starting encapsulated recursions..."

    perc = invariant_args["percentile"]
    distance_function = map_distance_function(invariant_args["metric"])

    data_dict = load_input_array_as_dictionary(representatives_path)
    n = length(data_dict)

    if !isnothing(cache_path) && isfile(cache_path)
        cache_dict = load_input_array_as_dictionary(cache_path)
    else
        cache_dict = copy(data_dict)
    end

    previous_groups = nothing
    if !isnothing(previous_groups_path) && isfile(previous_groups_path)
        previous_groups_dict = load_groups_as_dictionary(previous_groups_path)
        previous_groups = previous_groups_dict["groups"]
    end

    i = iteration
    full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
    while n > invariant_args["sample-size"]
        @info "Iteration: $i ($n objects)"
        flush(stdout)
        representatives_dict,cache_dict,groups = stratification_process(
            data_dict=data_dict,
            perc=perc,
            cache_dict=cache_dict,
            previous_groups=previous_groups,
            distance_function=distance_function,
            label=full_label,
            output_dir=invariant_args["output-dir"]
        )
        previous_groups = groups
        data_dict = representatives_dict
        n = length(data_dict)
        @info "\t$n objects after recursive iteration."
        i += 1
        full_label = "$(label).iteration_$(lpad(string(i),8,'0'))"
    end
    return
end

function stratification_process(;data_dict,perc,cache_dict,previous_groups,distance_function,label,output_dir)

    start_time = time()

    n_comparisons = 0
    D,index_map,threshold,x = compute_pairwise_distance_matrix(data_dict,perc,distance_function)
    n_comparisons += x

    representative_dict,groups = stratification_precomputed(data_dict,D,index_map,threshold)
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
    n_comparisons += x
    direct_groupsize_dict = Dict( id => length(group) for (id,group) in optimized_groups )

    if isnothing(previous_groups)
        previous_groups_label = "no"
    else
        previous_groups_label = "yes"
        optimized_groups = compile_previous_groupings(optimized_groups,previous_groups)
    end

    end_time = time()

    n = length(data_dict)
    G = length(optimized_dict)
    g = length(optimization_set)
    t = round(threshold,digits=4)
    N = sum(length(group) for group in values(optimized_groups))

    elapsed_time = round(end_time-start_time,digits=2)

    @info "\t[$label] $G groups ($g groups optimized); $n objects ($N total); threshold: $t (computed); $n_comparisons comparisons ($(nworkers())) CPU(s); cache: yes; previous_groups: $previous_groups_label); runtime: $(elapsed_time) seconds"

    groups_path = joinpath(output_dir,"$(label).groups.jsonl.gz")
    write_group_results(
        path=groups_path,
        label=label,
        stage="stratification_process",
        cache_label="yes",
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
    return optimized_dict,cache_dict,optimized_groups
end


