module GroupAggregation

    using Glob
    using Logging

    using ..PairwiseComparisons: map_distance_function,compute_pairwise_distance_matrix
    using ..FileHandling: load_milk_binaries,write_milk_binaries,load_groups_as_dictionary,
                           write_group_results,metadata_path,binary_path,group_meta_path,
                           cache_meta_path,cache_binary_path,count_rows,load_milk_cache
    using ..GroupStratification: ward_stratification,medoid_stratification,update_cache
    export group_aggregation

    function concatenate_files(;paths,concatenated_path)
        open(concatenated_path,"w") do io
            for path in paths
                write(io,read(path))
                rm(path)
            end
        end
    end

    function concatenating_aggregation(label,partition_dir,invariant_args)
        concat_representatives_path = joinpath(invariant_args["output-dir"],"$(label).concatenated.representatives.ids")
        representatives_paths = sort(glob("*.representatives.ids",partition_dir))
        concatenate_files(paths=representatives_paths,concatenated_path=concat_representatives_path)
        concatenate_files(paths=binary_path.(representatives_paths),concatenated_path=binary_path(concat_representatives_path))
        concatenate_files(paths=group_meta_path.(representatives_paths),concatenated_path=group_meta_path(concat_representatives_path))
        if isfile(cache_meta_path(first(representatives_paths)))
            concatenate_files(paths=cache_meta_path.(representatives_paths),concatenated_path=cache_meta_path(concat_representatives_path))
            concatenate_files(paths=cache_binary_path.(representatives_paths),concatenated_path=cache_binary_path(concat_representatives_path))
        end

        concat_groups_path = joinpath(invariant_args["output-dir"],"$(label).concatenated.groups.jsonl.gz")
        concatenate_files(paths=sort(glob("*.groups.jsonl.gz",partition_dir)),concatenated_path=concat_groups_path)
        concatenate_files(paths=sort(glob("*.groups.metadata.tsv",partition_dir)),concatenated_path=metadata_path(concat_groups_path))
        n_groups = sum(parse(Int,split(line,'\t')[6]) for line in eachline(metadata_path(concat_groups_path)); init=0)

        return concat_representatives_path,concat_groups_path,n_groups
    end

    """
    Aggregation of groups across partitions.
    """
    function group_aggregation(;label,threshold,input_path,partition_dir,invariant_args)
        concat_representatives_path,concat_groups_path,n_groups = concatenating_aggregation(label,partition_dir,invariant_args)
        @info "\t$n_groups groups after concatenation across partitioned results."
        flush(stdout)
        
        representatives_path = joinpath(invariant_args["output-dir"],"$(label).representatives.ids")
        groups_path = joinpath(invariant_args["output-dir"],"$(label).groups.jsonl.gz")
        if count_rows(concat_representatives_path) <= invariant_args["merge-threshold"]
            stratifying_aggregation(
                representatives_path=representatives_path,
                groups_path=groups_path,
                concat_representatives_path=concat_representatives_path,
                concat_groups_path=concat_groups_path,
                threshold=threshold,
                linkage=invariant_args["linkage"],
                percentile=invariant_args["percentile"],
                label=label,
                metric=invariant_args["metric"],
                cache_size=invariant_args["cache-size"],
                input_path=input_path
            )
        else
            mv(concat_representatives_path,representatives_path)
            mv(binary_path(concat_representatives_path),binary_path(representatives_path))
            mv(group_meta_path(concat_representatives_path),group_meta_path(representatives_path))
            if isfile(cache_meta_path(concat_representatives_path))
                mv(cache_meta_path(concat_representatives_path),cache_meta_path(representatives_path))
                mv(cache_binary_path(concat_representatives_path),cache_binary_path(representatives_path))
            end
            mv(concat_groups_path,groups_path)
            mv(metadata_path(concat_groups_path),metadata_path(groups_path))
        end
        return representatives_path,groups_path
    end

    function stratifying_aggregation(;representatives_path,groups_path,concat_representatives_path,concat_groups_path,threshold,
                                        linkage,percentile,label,metric,cache_size,input_path)
    
        start_time = time()
        distance_function = map_distance_function(metric)
        
        concat_representatives_dict,concat_groupsize_dict = load_milk_binaries(concat_representatives_path)
        cache_dict = cache_size >= 1 ? load_milk_cache(concat_representatives_path,concat_representatives_dict,concat_groupsize_dict) : nothing
        # partition representatives' own caches from the level input (the caches their radii were measured against),
        # so radius + merge distance bounds each member's distance to the final representative
        input_cache_dict = cache_size >= 1 ? load_milk_cache(input_path,concat_representatives_dict,concat_groupsize_dict;keep_ids=Set(keys(concat_representatives_dict))) : nothing
        concat_groups_info = load_groups_as_dictionary(concat_groups_path) # partition rep => members, radius
        concat_groups = concat_groups_info["groups"]

        D,index_map,_,n_comparisons = compute_pairwise_distance_matrix(concat_representatives_dict,percentile,distance_function,false;cache_dict=input_cache_dict)
        if startswith(linkage,"ward")
            balanced = linkage == "ward_balanced"
            optimized_dict,optimized_groups,distances_dict,_ = ward_stratification(concat_representatives_dict,concat_groupsize_dict,D,index_map,threshold,percentile,balanced)
            specificity_dict = Dict(id => Int32[] for id in keys(optimized_groups))
        else
            optimized_dict,optimized_groups,distances_dict,specificity_dict = medoid_stratification(concat_representatives_dict,concat_groupsize_dict,D,index_map,threshold;radius_dict=concat_groups_info["radius"])
        end
        
        total_groupsize_dict = Dict(id => sum(concat_groupsize_dict[m] for m in group) for (id,group) in optimized_groups) # members = partition reps
        cache_dict = cache_size >= 1 ? update_cache(optimized_dict,optimized_groups,cache_dict,cache_size,distance_function) : nothing
    
        # each previous-level member's distance to its merged representative, aligned with the expanded group:
        # dist(member, partitioned rep) + dist(partitioned rep, merged rep) is an upper bound (triangle inequality)
        partition_member_dists = concat_groups_info["distances"] # partitioned rep => its members' distances to it, in group order
        merged_distances_dict = Dict{String,Vector{Float32}}()
        for (merged_rep,partitioned_reps) in optimized_groups
            updated_dists = Float32[]
            for (i,partitioned_rep) in enumerate(partitioned_reps)
                if startswith(linkage,"ward")
                    # ward lists every member, representative included
                    offset = distances_dict[merged_rep][i] # partition centroid -> merged centroid
                    member_dists = partition_member_dists[partitioned_rep]
                    if isempty(member_dists)
                        # single-group partitions record no distances; the rep is its own centroid
                        member_dists = zeros(length(concat_groups[partitioned_rep]))
                    end
                    for member_dist in member_dists
                        push!(updated_dists,member_dist+offset)
                    end
                else
                    # medoid lists members excluding the representative
                    offset = D[index_map[partitioned_rep],index_map[merged_rep]]
                    j = 0 # position in partition_member_dists, which skips the partitioned rep
                    for member in concat_groups[partitioned_rep]
                        if member == partitioned_rep
                            if partitioned_rep != merged_rep
                                push!(updated_dists,offset)
                            end
                        else
                            j += 1
                            push!(updated_dists,partition_member_dists[partitioned_rep][j]+offset)
                        end
                    end
                end
            end
            merged_distances_dict[merged_rep] = updated_dists
        end

        optimized_groups = Dict(id => reduce(vcat,(concat_groups[m] for m in group)) for (id,group) in optimized_groups)
        direct_groupsize_dict = Dict(id => length(group) for (id,group) in optimized_groups) # members = previous level

        G = length(optimized_dict)
        N = sum(values(total_groupsize_dict))

        @info "\t[MERGE: $label] $G groups; $N objects; threshold: $(round(threshold,digits=4)); runtime: $(round(time()-start_time,digits=2)) seconds"
        flush(stdout)
    
        write_milk_binaries(optimized_dict,total_groupsize_dict,representatives_path;cache_dict=cache_dict)
        write_group_results(
            path=groups_path,
            label=label,
            stage="representative_stratification_merge",
            n_cached=count_rows(representatives_path),
            n_input_objects=N,
            n_groups=G,
            groups=optimized_groups,
            direct_groupsize_dict=direct_groupsize_dict,
            total_groupsize_dict=total_groupsize_dict,
            distances_dict=merged_distances_dict,
            specificity_dict=Dict(id => Int32[] for id in keys(optimized_groups)),
            threshold=threshold,
            n_comparisons=n_comparisons,
            cache_dict=cache_dict
        )
        rm(concat_representatives_path)
        rm(concat_groups_path)
        rm(metadata_path(concat_groups_path))
        rm(binary_path(concat_representatives_path))
        rm(group_meta_path(concat_representatives_path))
        rm(cache_meta_path(concat_representatives_path),force=true)
        rm(cache_binary_path(concat_representatives_path),force=true)
    end
end