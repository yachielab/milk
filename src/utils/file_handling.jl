using JSON
using Glob
using CodecZlib
using Statistics

function load_groups_as_dictionary(path)
    groups = Dict{String,Vector{String}}()
    spread_dict = Dict{String,Any}()
    specificity_dict = Dict{String,Any}()
    resolution_dict = Dict{String,Any}()
    thresholds = Vector{Float32}()
    open_file_read(path,gzip=true) do file
        for line in eachline(file)
            info = JSON.parse(line)
            representative_id = info["representative_id"]
            groups[representative_id] = info["group"]
            spread_dict[representative_id] = isempty(info["distances"]) ? nothing : mean(info["distances"]) 
            specificity_dict[representative_id] = isempty(info["specificity"]) ? nothing : mean(info["specificity"]) 
            resolution_dict[representative_id] = length(info["distances"])
            push!(thresholds,info["threshold"])
        end
    end
    groups_dict = Dict(
        "groups" => groups,
        "thresholds" => thresholds,
        "spread" => spread_dict,
        "specificity" => specificity_dict,
        "resolution" => resolution_dict
    )
    return groups_dict
end

function write_group_results(;path,label,stage,cache_label,compiled_label,n_input_objects,n_groups,groups,optimization_set,direct_groupsize_dict,distances_dict,specificity_dict,threshold,n_comparisons)
    """
    """
    open_file_write(path,gzip=true) do file
        for (representative_id,group) in groups
            group_info = Dict(
                "label" => label,
                "stage" => stage,
                "cache" => cache_label,
                "compiled" => compiled_label,
                "representative_id" => representative_id,
                "group" => group,
                "direct_group_size" => direct_groupsize_dict[representative_id],
                "compiled_group_size" => length(group),
                "distances" => distances_dict[representative_id],
                "specificity" => specificity_dict[representative_id],
                "optimized" => (representative_id in optimization_set),
                "threshold" => threshold,
                "n_input_objects" => n_input_objects,
                "n_groups" => n_groups,
                "n_comparisons" => n_comparisons
            )
            JSON.print(file,group_info)
            println(file)
        end
    end
    return
end

function write_dictionary_as_csv(dict,path,gzip=false)
    open_file_write(path,gzip=gzip) do handle
        for (id,vec) in dict
            write(handle,"$id,$(join(vec,','))\n")
        end
    end
end

function load_values_as_list(path,dtype)
    if dtype == "float"
        vals = Vector{Float64}()
        open(path,"r") do handle
            for val in readlines(handle)
                push!(vals,parse(Float32,val))
            end
        end
    elseif dtype == "integer"
        vals = Vector{Int64}()
        open(path,"r") do handle
            for val in readlines(handle)
                push!(vals,parse(Int,val))
            end
        end   
    elseif dtype == "string"
        vals = Vector{String}()
        open(path,"r") do handle
            for val in readlines(handle)
                push!(vals,val)
            end
        end
    else
        throw("Unrecognized data type specified!")
    end
    return vals
end

function load_input_array_as_lists(path)
    ids = Vector{String}()
    vecs = Vector{Vector{Float32}}()
    open_file_read(path,gzip=false) do handle
        for line in readlines(handle)
            entry = split(line,",")
            id = entry[1]
            # vec = parse.(Float32,entry[2:end])
            vec = [val == "" ? NaN32 : parse(Float32,val) for val in entry[2:end]]
            if !any(isnan,vec)
                push!(ids,id)
                push!(vecs,vec)
            end
        end
    end
    return ids,vecs
end

function load_input_array_as_dictionary(path)
    input_dict = Dict{String,Vector{Float32}}()
    open_file_read(path,gzip=false) do handle
        for line in readlines(handle)
            entry = split(line,",")
            id = entry[1]
            # vec = parse.(Float32,entry[2:end])
            vec = [val == "" ? NaN32 : parse(Float32,val) for val in entry[2:end]]
            if !any(isnan,vec)
                input_dict[id] = vec
            end
        end
    end
    return input_dict
end

function partition_input_file(input_path,label,invariant_args)

    partition_dir = joinpath(invariant_args["output-dir"],"$(label).split")
    mkpath(partition_dir) # does not throw error if it exists

    function get_partitioned_input_path(p)
        partition_label = "partition_$(lpad(string(p),8,'0'))"
        return joinpath(partition_dir, "$(label).$(partition_label).csv")
    end

    p = 1
    buffer = []
    open_file_read(input_path,gzip=false) do instream
        for line in eachline(instream)
            if length(buffer) == invariant_args["partition-size"]
                partitioned_input_path = get_partitioned_input_path(p)
                open_file_write(partitioned_input_path,gzip=false) do outstream
                    write(outstream,join(buffer,"\n")*"\n")
                end
                empty!(buffer)
                p += 1
            end
            push!(buffer,line)
        end
    end

    if !isempty(buffer)
        partitioned_input_path = get_partitioned_input_path(p)
        open_file_write(partitioned_input_path,gzip=false) do outstream
            write(outstream,join(buffer,"\n")*"\n")
        end
    end
    return partition_dir
end

function batch_partitioned_files(partition_dir,label,invariant_args)

    function get_batch_directory(b)
        batch_label = "batch_$(lpad(string(b),8,'0'))"
        return joinpath(partition_dir,"$(label).$(batch_label).work")
    end

    pattern = "*.csv"
    paths = sort(glob(pattern,partition_dir))

    files = []
    batches = []
    if length(paths) > invariant_args["batch-size"]
        for (b,batch) in enumerate(Iterators.partition(paths,invariant_args["batch-size"]))
            batch_dir = get_batch_directory(b)
            mkpath(batch_dir)
            push!(batches,batch_dir)
            for path in batch
                updated_path = joinpath(batch_dir,basename(path))
                push!(files,updated_path)
                mv(path,updated_path)
            end
        end
    else
        batch_dir = get_batch_directory(0)
        mkpath(batch_dir)
        push!(batches,batch_dir)
        for path in paths
            updated_path = joinpath(batch_dir,basename(path))
            push!(files,updated_path)
            mv(path,updated_path)
        end
    end
    return files,batches
end

function partition_and_batch_input_files(input_path,label,invariant_args)
    partition_dir = partition_input_file(input_path,label,invariant_args)
    input_paths,batch_dirs = batch_partitioned_files(partition_dir,label,invariant_args)
    return partition_dir,input_paths,batch_dirs
end

function get_object_count(path)
    command = pipeline(`cat $path`,`wc -l`)
    return parse(Int,strip(read(command,String)))
end

function open_file_write(f::Function, path; gzip=true)
    stream = gzip ? GzipCompressorStream(open(path,"w")) : open(path,"w")
    try
        return f(stream)  # Pass the stream to the function
    finally
        close(stream)  # Ensure the stream is closed properly
    end
end

function open_file_read(f::Function, path; gzip=true)
    stream = gzip ? GzipDecompressorStream(open(path,"r")) : open(path, "r")
    try
        return f(stream)  # Pass the stream to the function
    finally
        close(stream)  # Ensure the stream is closed properly
    end
end

function concatenate_files(paths,concatenated_path; gzip=false)
    x = 0
    open_file_write(concatenated_path,gzip=gzip) do outstream
        for path in paths
            open_file_read(path,gzip=gzip) do instream
                for line in eachline(instream)
                    println(outstream,line)
                    x += 1
                end
            end
            rm(path)
        end
    end
    return x
end

function concatenate_partitioned_results(label,partition_dir,invariant_args)
    representatives_pathlist = glob("*.representatives.csv",partition_dir)
    groups_pathlist = glob("*.groups.jsonl.gz",partition_dir)
    concat_representatives_path = joinpath(invariant_args["output-dir"],"$(label).concatenated.representatives.csv")
    concat_groups_path = joinpath(invariant_args["output-dir"],"$(label).concatenated.groups.jsonl.gz")
    concatenate_files(representatives_pathlist,concat_representatives_path,gzip=false)
    x = concatenate_files(groups_pathlist,concat_groups_path,gzip=true)
    return concat_representatives_path,concat_groups_path,x
end

function is_broken_symlink(path)
    return islink(path) && !isfile(path)
end

function clean_directory(work_dir,partition_dir,exclusion_set)
    rm(partition_dir,recursive=true)
    pattern = "*.representatives.csv"
    for path in glob(pattern,work_dir)
        if path in exclusion_set continue end
        rm(path)
    end
    pattern = "*.input.csv"
    for path in glob(pattern,work_dir)
        if is_broken_symlink(path)
            rm(path)
        end
    end
    return
end

function write_value_as_txt(value,path)
    open(path,"w") do handle
        write(handle,"$value\n")
    end
end

function write_values_as_txt(values,path)
    open(path,"w") do handle
        for value in values
            write(handle,"$value\n")
        end
    end
end
