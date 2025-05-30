using GZip
using Glob
using CodecZlib
using Distances
using Base.Iterators
using Base.Threads
using Logging

include("file_handling.jl")

function group_stratification_precomputed_distances(;ids,vecs,D,index_map,threshold)
    representative_dict = Dict{String,Vector{Float32}}()
    groups = Dict{String,Vector{String}}()
    for (candidate_id,candidate_vec) in zip(ids,vecs)
        i = index_map[candidate_id]
        closest_id   = nothing
        closest_dist = Inf
        for (representative_id,representative_vec) in representative_dict
            j = index_map[representative_id]
            @inbounds dist = D[i,j]
            if (dist <= threshold) && (dist < closest_dist)
                closest_id   = representative_id
                closest_dist = dist
            end
        end
        if isnothing(closest_id)
            representative_dict[candidate_id] = candidate_vec
            groups[candidate_id] = [candidate_id]
        else
            push!(groups[closest_id],candidate_id)
        end
    end
    # @info "\tgroup_stratification_precomputed: $(length(groups))..."
    return representative_dict,groups
end

function group_stratification_predefined_medoids(;representative_dict,ids,vecs,D,index_map,cache_dict,threshold,distance_function)
    groups = Dict( id => [id] for id in keys(representative_dict) )
    distances_dict = Dict{String,Float16}()
    specificity_dict = Dict{String,Int32}()
    n_comparisons = 0
    for (candidate_id,candidate_vec) in zip(ids,vecs)
        if haskey(representative_dict,candidate_id) continue end
        i = index_map[candidate_id]
        closest_id = nothing
        closest_dist = Inf
        specificity = 0
        for (representative_id,representative_vec) in representative_dict
            if haskey(index_map,representative_id)
                j = index_map[representative_id]
                @inbounds dist = D[i,j]
            elseif !isnothing(cache_dict) && (haskey(cache_dict,representative_id))
                dist = distance_function(candidate_vec,representative_vec)
                n_comparisons += 1
            else
                # should never hit this
                error("Representative id ($representative_id) found in neither the index_map nor cache_dict...")
            end
            specificity += Int(dist <= threshold)
            if dist < closest_dist
                closest_id = representative_id
                closest_dist = dist
            end
        end
        push!(groups[closest_id],candidate_id)
        distances_dict[candidate_id] = closest_dist
        specificity_dict[candidate_id] = specificity
    end
    return groups,distances_dict,specificity_dict,n_comparisons
end

function compute_centroid(medoid_dict)
    arr = reduce(hcat,values(medoid_dict))
    return vec(mean(arr,dims=2))
end

function representative_optimization(;ids,vecs,index_map,groups,optimization_set,cache_dict,distance_function)
    optimized_dict= Dict{String,Vector{Float32}}()
    n_comps = 0
    for (representative_id,group) in groups
        if representative_id in optimization_set
            if (length(group) == 1) # should never encounter if optimization_set is set properly -- remove after test
                error("Unexpected state: group ($representative_id) has a size of 1 -- unable to optimize!")
            end
            medoid_dict  = get_medoid_dictionary_from_vectors(group,ids,vecs,index_map,cache_dict)
            centroid_vec = compute_centroid(medoid_dict)
            medoid_id,medoid_vec,n_centroid_comps = determine_medoid(medoid_dict,centroid_vec,distance_function)
            n_comps += n_centroid_comps
            optimized_dict[medoid_id] = medoid_vec
        else
            optimized_dict[representative_id] = vecs[index_map[representative_id]]
        end
    end
    return optimized_dict,n_comps
end

function get_medoid_dictionary_from_vectors(group,ids,vecs,index_map,cache_dict)
    medoid_dict = Dict{String,Vector{Float32}}() # this is for a single group
    for sample_id in group
        if haskey(index_map,sample_id)
            medoid_dict[sample_id] = vecs[index_map[sample_id]]
        elseif !isnothing(cache_dict) && haskey(cache_dict,sample_id)
            medoid_dict[sample_id] = cache_dict[sample_id]
        end
    end
    return medoid_dict
end

function get_medoid_dictionary_from_representatives(group,representative_dict,cache_dict)
    medoid_dict = Dict{String,Vector{Float32}}() # this is for a single group
    for sample_id in group
        if haskey(representative_dict,sample_id)
            medoid_dict[sample_id] = representative_dict[sample_id]
        elseif !isnothing(cache_dict) && haskey(cache_dict,sample_id)
            medoid_dict[sample_id] = cache_dict[sample_id]
        end
    end
    return medoid_dict
end

function determine_medoid(medoid_dict,centroid_vec,distance_function)
    medoid_id = nothing
    medoid_vec = nothing
    medoid_dist = Inf
    n_comparisons = 0
    for (sample_id,vec) in medoid_dict
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

function optimize_compiled_representatives(representative_dict,groups,cache_dict,metric,optimization_list)
    optimization_set = Set(optimization_list)
    optimized_dict   = Dict{String,Vector{Float32}}()
    optimized_groups = Dict{String,Vector{String}}()
    n = 0
    for (representative_id,group) in groups
        if representative_id in optimization_set
            medoid_dict = get_medoid_dict(group,representative_dict,cache_dict)
            centroid_vec = vec(mean(reduce(hcat,values(medoid_dict)),dims=2))
            medoid_id,medoid_vec,m = determine_medoid(medoid_dict,centroid_vec,metric)
            optimized_dict[medoid_id] = medoid_vec
            optimized_groups[medoid_id] = group
            n += m
        else
            optimized_dict[representative_id] = representative_dict[representative_id]
            optimized_groups[representative_id] = group
        end
    end
    return optimized_dict,optimized_groups,n
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
    optimization_set = Set{String}()
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
            push!(optimization_set,representative_id)
            compiled_groups[representative_id] = Vector{String}()
            for id in group
                if !haskey(previous_groups, id)
                    error("Missing previous group for id: $id")
                end
                append!(compiled_groups[representative_id],previous_groups[id])
            end
        end
    end
    return compiled_groups,optimization_set
end

function merge_partitioned_results(concat_files,label,cache_path,previous_groups_path,invariant_args)
    concat_n_groups_path        = concat_files[1]
    concat_thresholds_path      = concat_files[2]
    concat_n_comparisons_path   = concat_files[3]
    concat_representatives_path = concat_files[4]
    concat_groups_paths         = concat_files[5]

    m = sum(load_values_to_list(concat_n_groups_path,"integer"))

    n_groups_path        = joinpath(invariant_args["output-dir"],"$(label).n_groups.txt")
    thresholds_path      = joinpath(invariant_args["output-dir"],"$(label).thresholds.txt")
    n_comparisons_path   = joinpath(invariant_args["output-dir"],"$(label).n_comparisons.txt")
    representatives_path = joinpath(invariant_args["output-dir"],"$(label).representatives.csv.gz")
    groups_path          = joinpath(invariant_args["output-dir"],"$(label).groups.jsonl.gz")
    if m <= MERGE_THRESHOLD
        command_list = [
            "julia",
            "merge_partitioned_results.jl",
            "-i", concat_representatives_path,
            "-g", concat_groups_paths,
            "-t", concat_thresholds_path,
            "-n", concat_n_comparisons_path,
            "-l", label,
            "-c", cache_path,
            "-m", invariant_args["metric"],
            "-r", invariant_args["seed"],
            "-o", invariant_args["output-dir"]
        ]
        if !isnothing(previous_groups_path)
            push!(command_list,"-G")
            push!(command_list,previous_groups_path)
        end
        command = Cmd(command_list)
        run(command)
    else
        mv(concat_n_groups_path,n_groups_path)
        mv(concat_thresholds_path,thresholds_path)
        mv(concat_n_comparisons_path,n_comparisons_path)
        mv(concat_representatives_path,representatives_path)
        mv(concat_groups_path,groups_path)
    end
    return (n_groups_path,thresholds_path,n_comparisons_path,representatives_path,groups_path)
end

function stratify_partitioned_inputs(;batches,files,label,cache_path,previous_groups_path,partition_dir,invariant_args)
    b = length(batches)
    if b == 1
        # @info "\tSingle batch detected: direct execution..."
        for partitioned_input_path in files
            partition_label = replace(basename(partitioned_input_path), ".csv.gz" => "")
            command = [
                "julia",
                "--threads", string(invariant_args["threads"]),
                "threshold_based_group_stratification.jl",
                "-i", partitioned_input_path,
                "-t", string(invariant_args["percentile"]),
                "-l", partition_label,
                "-m", invariant_args["metric"],
                "-o", partition_dir
            ]
            if !isnothing(previous_groups_path)
                push!(command,"-p",previous_groups_path)
            end
            if !isnothing(cache_path)
                push!(command,"-c",cache_path)
            end
            command = Cmd(command)
            run(command)
        end
    else
        println("\tSubmitting $b batches as an array job...")
        flush(stdout)
        submit_stratification_batch_array_jobs(
            batches=batches,
            files=files,
            label=label,
            partition_dir=partition_dir,
            percentile=invariant_args["percentile"],
            metric=invariant_args["metric"],
            threads=invariant_args["threads"]
        )
    end
    return
end

