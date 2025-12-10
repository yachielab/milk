module HierarchicalReconstruction

    using Glob
    using JSON
    using Logging
    using Statistics
    using ..FileHandling: open_file_write,open_file_read

    export hierarchical_reconstruction

    function extract_iteration(path)
        iteration = replace(basename(path), ".gz" => "", ".groups.jsonl" => "")
        iteration = split(split(iteration,".")[end],"_")[end]
        return parse(Int,iteration)
    end

    @inline function get_spread(group_dict)
        """ group spread (inverse density) of members with respect to representative
        """
        return isempty(group_dict["distances"]) ? "" : string(mean(group_dict["distances"]))
    end

    @inline function get_specificity(group_dict)
        """ group specificity of members with respect to representative
        """
        return isempty(group_dict["specificity"]) ? "" : string(mean(group_dict["specificity"]))
    end

    function hierarchical_reconstruction(input_dir)

        @info "Hierarchical reconstruction"
        flush(stdout)

        output_dir = joinpath(input_dir,"output")
        mkpath(output_dir) # exist_ok
        vertices_path = joinpath(output_dir,"vertices.csv.gz")
        edges_path = joinpath(output_dir,"edges.csv.gz")

        pathlist = sort(glob("*.groups.jsonl.gz",input_dir))
        if length(pathlist) == 0
            error("No '.groups.jsonl.gz' files found in directory: $(input_dir)")
        end

        id_map = Dict{String,String}()
        id = 1
        open_file_write(vertices_path, gzip=true) do vertices_io
            open_file_write(edges_path, gzip=true) do edges_io
                println(vertices_io,"node_id,representative_id,group_size,iteration,threshold,spread,specificity,resolution")
                println(edges_io,"source,target")

                groups_path = pathlist[1]
                i = extract_iteration(groups_path)+1
                open_file_read(groups_path,gzip=true) do file
                    for line in eachline(file)
                        group_dict = JSON.parse(line)
                        group_label = string(group_dict["representative_id"],"|",i)
                        id_map[group_label] = "I$(id)"
                        id += 1
                        group_id = id_map[group_label]

                        # instantiating leaves
                        for sample_id in group_dict["group"]
                            println(vertices_io,"$(sample_id),$(sample_id),1,0,0,,,0")
                            println(edges_io,"$(group_id),$(sample_id)")
                        end
        
                        fields = [
                            group_id,
                            group_dict["representative_id"],
                            group_dict["compiled_group_size"],
                            i,
                            group_dict["threshold"],
                            get_spread(group_dict),
                            get_specificity(group_dict),
                            length(group_dict["distances"]) # resolution
                        ]
                        println(vertices_io,join(fields,","))
                    end
                end

                for groups_path in pathlist[2:end]
                    i = extract_iteration(groups_path)+1
                    open_file_read(groups_path,gzip=true) do file
                        for line in eachline(file)
                            group_dict = JSON.parse(line)
                            group_label = string(group_dict["representative_id"],"|",i)
                            id_map[group_label] = "I$(id)"
                            id += 1
                            group_id = id_map[group_label]

                            fields = [
                                group_id,
                                group_dict["representative_id"],
                                group_dict["compiled_group_size"],
                                i,
                                group_dict["threshold"],
                                get_spread(group_dict),
                                get_specificity(group_dict),
                                length(group_dict["distances"]) # resolution
                            ]
                            println(vertices_io,join(fields,","))
                            for sample_id in group_dict["group"]
                                subgroup_label = "$(sample_id)|$(i-1)"
                                if !haskey(id_map,subgroup_label) continue end
                                subgroup_id = id_map[subgroup_label]
                                println(edges_io,"$(group_id),$(subgroup_id)")
                            end
                        end
                    end
                end
            end
        end
    end

end