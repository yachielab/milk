using GZip
using Glob
using CodecZlib
using Distances
using Base.Iterators
using Base.Threads
using Logging

include("hpc.jl")
include("file_handling.jl")

function map_object_precomputed(i,representative_dict,D,index_map,threshold)
    closest_id = nothing
    closest_dist = Inf
    specificity = 0
    for (representative_id,representative_vec) in representative_dict
        j = index_map[representative_id]
        @inbounds dist = D[i,j]
        if dist <= threshold
            specificity += 1
            if dist < closest_dist
                closest_id = representative_id
                closest_dist = dist
            end
        end
    end
    return closest_id,closest_dist,specificity
end

function stratification_precomputed(data_dict,D,index_map,threshold)
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

function map_object(vec,representative_dict,threshold,distance_function)
    closest_id = nothing
    closest_dist = Inf
    specificity = 0
    n_comparisons = 0
    for (representative_id,representative_vec) in representative_dict
        dist = distance_function(vec,representative_vec)
        n_comparisons += 1
        if dist <= threshold
            specificity += 1
            if dist < closest_dist
                closest_id = representative_id
                closest_dist = dist
            end
        end
    end
    return closest_id,closest_dist,n_comparisons,specificity
end

function stratification_predefined_medoids_precomputed(;data_dict,representative_dict,D,index_map,cache_dict,threshold,distance_function)
    groups = Dict( id => [id] for id in keys(representative_dict) )
    distances_dict = Dict(id => Vector{Float16}() for id in keys(representative_dict))
    specificity_dict = Dict(id => Vector{Int32}() for id in keys(representative_dict))
    for (candidate_id,candidate_vec) in data_dict
        if haskey(representative_dict,candidate_id) continue end
        i = index_map[candidate_id]
        closest_id,closest_dist,specificity = map_object_precomputed(i,representative_dict,D,index_map,threshold)
        if isnothing(closest_id)
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

function compute_centroid(medoid_dict)
    return vec(mean(reduce(hcat,values(medoid_dict)),dims=2))
end

function optimize_representatives(data_dict,groups,cache_dict,distance_function)
    optimized_dict = Dict{String,Vector{Float32}}()
    optimization_set = Set{String}()
    x = 0
    for (representative_id,group) in groups
        if length(group) == 1
            optimized_dict[representative_id] = data_dict[representative_id]
        else
            medoid_dict,cache_set = get_medoid_dictionary(group,data_dict,cache_dict)
            if length(medoid_dict) < 3
                optimized_dict[representative_id] = data_dict[representative_id]
            else
                centroid_vec = compute_centroid(medoid_dict)
                medoid_id,medoid_vec,n_centroid_comps = determine_medoid(medoid_dict,centroid_vec,cache_set,distance_function)
                x += n_centroid_comps
                optimized_dict[medoid_id] = medoid_vec
                push!(optimization_set,medoid_id)
            end
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

function determine_medoid(medoid_dict,centroid_vec,cache_set,distance_function)
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


function attempt_to_merge_partitioned_results(;concat_representatives_path,concat_groups_path,p,m,label,cache_path,invariant_args)
    representatives_path = joinpath(invariant_args["output-dir"],"$(label).representatives.csv.gz")
    groups_path = joinpath(invariant_args["output-dir"],"$(label).groups.jsonl.gz")
    if (p > 1) && (m <= invariant_args["merge-threshold"])
        # @info "\tMerging: number of concatenated groups, $m <= merge threshold ($(invariant_args["merge-threshold"]))"
        command = [
            "julia",
            joinpath(dirname(@__DIR__),"representative_stratification_merge.jl"),
            "-i", concat_representatives_path,
            "-g", concat_groups_path,
            "-l", label,
            "-m", invariant_args["metric"],
            "--verbose",
            "-o", invariant_args["output-dir"]
        ]
        if !isnothing(cache_path)
            push!(command,"-c",cache_path)
        end
        command = Cmd(command)
        run(command)
        # merged_representatives_path = joinpath(invariant_args["output-dir"],"$(label).merged.representatives.csv.gz")
        # merged_groups_path = joinpath(invariant_args["output-dir"],"$(label).merged.groups.jsonl.gz")
        # symlink(merged_representatives_path,representatives_path)
        # symlink(merged_groups_path,groups_path)
        # TODO: keep merged and concatenated in work_dir and just mv final groups file into out_dir (no more symlink)
    else
        mv(concat_representatives_path,representatives_path)
        mv(concat_groups_path,groups_path)
    end
    return representatives_path,groups_path
end

function stratify_partitioned_inputs(;batches,files,label,cache_path,previous_groups_path,partition_dir,invariant_args)
    b = length(batches)
    if b == 1
        # @info "\tSingle batch detected: direct execution..."
        for partitioned_input_path in files
            partition_label = split(replace(basename(partitioned_input_path),".csv.gz" => ""),".")[end]
            command = [
                "julia",
                "--threads", string(invariant_args["threads"]),
                joinpath(dirname(@__DIR__),"threshold_based_group_stratification.jl"),
                "-i", partitioned_input_path,
                "-t", string(invariant_args["percentile"]),
                "-l", partition_label,
                "-m", invariant_args["metric"],
                "--verbose",
                "-o", partition_dir
            ]
            if !isnothing(cache_path)
                push!(command,"-c",cache_path)
            end
            if !isnothing(previous_groups_path)
                push!(command,"-G",previous_groups_path)
            end
            command = Cmd(command)
            run(command)
        end
    else
        @info "\tSubmitting $b batches as an array job..."
        submit_stratification_batch_array_jobs(
            batches=batches,
            files=files,
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
    @info "\tInput file partitioned into $(length(files)) files ($(length(batches)) batches)"
    flush(stdout)

    # @info "\tstratifying partitioned input files..."
    stratify_partitioned_inputs(
        batches=batches,
        files=files,
        label=label,
        cache_path=cache_path,
        previous_groups_path=previous_groups_path,
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
        label=label,
        cache_path=cache_path,
        invariant_args=invariant_args
    )
    delete_work_directory(partition_dir)
    return representatives_path,groups_path
end
