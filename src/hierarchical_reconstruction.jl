using ArgParse
using Logging
using Glob
using JSON

include("utils/file_handling.jl")

function parse_arguments()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--input-dir","-i"
            arg_type = String
            required = true
            help = ""
    end
    return parse_args(args)
end

function extract_iteration(path)
    iteration = replace(basename(path), ".gz" => "", ".groups.jsonl" => "")
    iteration = split(split(iteration,".")[end],"_")[end]
    return parse(Int,iteration)
end

function main()
    args = parse_arguments()

    @info "Hierarchical reconstruction"
    flush(stdout)

    output_dir = joinpath(args["input-dir"],"output")
    try
        mkdir(output_dir)
    catch e
        @warn "Directory already exists: $e"
        flush(stdout)
    end
    vertices_path = joinpath(output_dir,"vertices.csv.gz")
    edges_path = joinpath(output_dir,"edges.csv.gz")

    pathlist = sort(glob("*.groups.jsonl.gz",args["input-dir"]))
    if length(pathlist) == 0
        error("No '.groups.jsonl.gz' files found in directory: $(args["input-dir"])")
    end

    id_map = Dict{String,String}()
    id = 1
    open_file_write(vertices_path, gzip=true) do vertices_io
        open_file_write(edges_path, gzip=true) do edges_io
            println(vertices_io,"node_id,representative_id,group_size,iteration,threshold,spread,specificity,resolution")
            println(edges_io,"source,target")

            groups_path = pathlist[1]
            iteration = extract_iteration(groups_path)+1
            open_file_read(groups_path,gzip=true) do file
                for line in eachline(file)
                    group_dict = JSON.parse(line)
                    group_label = string(group_dict["representative_id"],"|",iteration)
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
                        iteration,
                        group_dict["threshold"],
                        isempty(group_dict["distances"]) ? "" : string(mean(group_dict["distances"])), # spread (inverse density)
                        isempty(group_dict["specificity"]) ? "" : string(mean(group_dict["specificity"])), # specificity
                        length(group_dict["distances"]) # resolution
                    ]
                    println(vertices_io,join(fields,","))
                end
            end

            for groups_path in pathlist[2:end]
                iteration = extract_iteration(groups_path)+1
                open_file_read(groups_path,gzip=true) do file
                    for line in eachline(file)
                        group_dict = JSON.parse(line)
                        group_label = string(group_dict["representative_id"],"|",iteration)
                        id_map[group_label] = "I$(id)"
                        id += 1
                        group_id = id_map[group_label]

                        fields = [
                            group_id,
                            group_dict["representative_id"],
                            group_dict["compiled_group_size"],
                            iteration,
                            group_dict["threshold"],
                            isempty(group_dict["distances"]) ? "" : string(mean(group_dict["distances"])), # spread (inverse density)
                            isempty(group_dict["specificity"]) ? "" : string(mean(group_dict["specificity"])), # specificity
                            length(group_dict["distances"]) # resolution
                        ]
                        println(vertices_io,join(fields,","))
                        for sample_id in group_dict["group"]
                            subgroup_label = "$(sample_id)|$(iteration-1)"
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

main()
