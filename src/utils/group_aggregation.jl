module GroupAggregation

    using Glob
    using Logging

    using ..PairwiseComparisons: map_distance_function,map_object_to_representative
    using ..FileHandling: load_values_as_list,load_input_array_as_dictionary,load_groups_as_dictionary,
                           write_dictionary_as_csv,write_group_results,open_file_write,open_file_read,
                           attempt_to_load_cache
    using ..GroupStratification: stratification,stratification_predefined_medoids,compile_previous_groupings
    using ..RepresentativeOptimization: optimize_representatives

    export group_aggregation

    function concatenate_files(;paths,concatenated_path,gzip=false)
        n_groups = 0
        open_file_write(concatenated_path,gzip=gzip) do outstream
            for path in paths
                open_file_read(path,gzip=gzip) do instream
                    for line in eachline(instream)
                        println(outstream,line)
                        n_groups += 1
                    end
                end
                rm(path)
            end
        end
        return n_groups
    end

    function concatenating_aggregation(label,partition_dir,invariant_args)
        concat_representatives_path = joinpath(invariant_args["output-dir"],"$(label).concatenated.representatives.csv")
        concatenate_files(
            paths=sort(glob("*.representatives.csv",partition_dir)),
            concatenated_path=concat_representatives_path,
            gzip=false
        )
        concat_groups_path = joinpath(invariant_args["output-dir"],"$(label).concatenated.groups.jsonl.gz")
        n_groups = concatenate_files(
            paths=sort(glob("*.groups.jsonl.gz",partition_dir)),
            concatenated_path=concat_groups_path,
            gzip=true
        )
        return concat_representatives_path,concat_groups_path,n_groups
    end

    function group_aggregation(;label,threshold,cache_path,partition_dir,invariant_args)
        """
        Aggregation of groups across partitions.
        """
        concat_representatives_path,concat_groups_path,n_groups = concatenating_aggregation(label,partition_dir,invariant_args)
        @info "\t$n_groups groups after concatenation across partitioned results."
        flush(stdout)
        
        representatives_path = joinpath(invariant_args["output-dir"],"$(label).representatives.csv")
        groups_path = joinpath(invariant_args["output-dir"],"$(label).groups.jsonl.gz")
        if n_groups <= invariant_args["merge-threshold"]
            stratifying_aggregation(
                concat_representatives_path=concat_representatives_path,
                concat_groups_path=concat_groups_path,
                threshold=threshold,
                cache_path=cache_path,
                label=label,
                metric=invariant_args["metric"],
                output_dir=invariant_args["output-dir"]
            )
        else
            mv(concat_representatives_path,representatives_path)
            mv(concat_groups_path,groups_path)
        end
        return representatives_path,groups_path
    end

    function stratifying_aggregation(;concat_representatives_path,concat_groups_path,threshold,
                                                 cache_path,label,metric,output_dir)
        start_time = time()
        distance_function = map_distance_function(metric)

        concat_representatives_dict = load_input_array_as_dictionary(concat_representatives_path)
        concat_groups_dict = load_groups_as_dictionary(concat_groups_path)

        concat_groups = concat_groups_dict["groups"]
        n_comparisons = 0

        n = sum(length(group) for group in values(concat_groups))

        cache_dict = attempt_to_load_cache(cache_path)
        cache_label = isnothing(cache_dict) ? "no" : "yes"

        """ Additional round of group stratification over the representatives
        """
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

        @info "\t[MERGE: $label] $G groups ($g groups optimized); $n objects; threshold: $t; cache: $cache_label; runtime: $(round(end_time-start_time,digits=2)) seconds"
        flush(stdout)

        representatives_path = joinpath(output_dir,"$label.representatives.csv")
        groups_path = joinpath(output_dir,"$label.groups.jsonl.gz")
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
        rm(concat_representatives_path)
        rm(concat_groups_path)
    end
end