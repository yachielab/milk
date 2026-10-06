module FileHandling

    using JSON
    using Glob
    using Logging
    using CodecZlib
    using Statistics

    export load_groups_as_dictionary,
        write_group_results,
        partition_input_file,
        batch_partitioned_files,
        partition_and_batch_input_files,
        open_file_write,
        open_file_read,
        clean_directory,
        final_cleanup,
        write_values_as_txt,
        attempt_to_cache_file,
        metadata_path,
        binary_path,
        groupsize_path,
        stream_binary_file,
        load_milk_binaries,
        write_milk_binaries,
        convert_input_csv_to_binaries,
        prepare_milk_input


    ########## inline helpers ##########
    function metadata_path(groups_path)
        return replace(groups_path, r"\.jsonl\.gz$" => ".metadata.tsv")
    end

    function binary_path(ids_path)
        if !endswith(ids_path,".ids")
            error("Expected a .ids path: $(ids_path)")
        end
        return replace(ids_path,r"\.ids$" => ".bin")
    end

    function groupsize_path(ids_path)
        if !endswith(ids_path,".ids")
            error("Expected a .ids path: $(ids_path)")
        end
        return replace(ids_path,r"\.ids$" => ".sizes")
    end
    ####################################

    function load_groups_as_dictionary(path)

        threshold = split(readline(metadata_path(path)),'\t')[5]

        groups = Dict{String,Vector{String}}()
        spread_dict = Dict{String,Any}()
        specificity_dict = Dict{String,Any}()
        resolution_dict = Dict{String,Any}()
        open_file_read(path,gzip=true) do file
            for line in eachline(file)
                info = JSON.parse(line)
                representative_id = info["representative_id"]
                groups[representative_id] = info["group"]
                spread_dict[representative_id] = isempty(info["distances"]) ? nothing : mean(info["distances"]) 
                specificity_dict[representative_id] = isempty(info["specificity"]) ? nothing : mean(info["specificity"]) 
                resolution_dict[representative_id] = length(info["distances"])
            end
        end
        groups_dict = Dict(
            "groups" => groups,
            "threshold" => threshold,
            "spread" => spread_dict,
            "specificity" => specificity_dict,
            "resolution" => resolution_dict
        )
        return groups_dict
    end

    function write_group_metadata(;path,label,stage,cache_label,compiled_label,threshold,n_input_objects,n_groups,n_comparisons)
        open_file_write(metadata_path(path),gzip=false) do file
            # columns: label, stage, cache, compiled, threshold, n_input_objects, n_groups, n_comparisons
            fields = [label,stage,cache_label,compiled_label,threshold,n_input_objects,n_groups,n_comparisons]
            println(file,join(fields,'\t'))
        end
    end

    function write_group_results(;path,label,stage,cache_label,compiled_label,n_input_objects,n_groups,groups,optimization_set,direct_groupsize_dict,distances_dict,specificity_dict,threshold,n_comparisons)

        write_group_metadata(
            path=path,
            label=label,
            stage=stage,
            cache_label=cache_label,
            compiled_label=compiled_label,
            threshold=threshold,
            n_input_objects=n_input_objects,
            n_groups=n_groups,
            n_comparisons=n_comparisons
        )

        open_file_write(path,gzip=true) do file
            for (representative_id,group) in groups
                group_info = Dict(
                    "label" => label,
                    "representative_id" => representative_id,
                    "group" => group,
                    "direct_group_size" => direct_groupsize_dict[representative_id],
                    "total_group_size" => length(group),
                    "distances" => distances_dict[representative_id],
                    "specificity" => specificity_dict[representative_id],
                    "optimized" => (representative_id in optimization_set)
                )
                JSON.print(file,group_info)
                println(file)
            end
        end
        return
    end

    function stream_binary_file(f,ids_path)
        n = countlines(ids_path)
        if n == 0
            return
        end
        d = div(filesize(binary_path(ids_path)),4*n) # vector dimensionality (number of Float32 per row)
        vec = Vector{Float32}(undef,d)
        open(binary_path(ids_path)) do bin_io
            for id in eachline(ids_path)
                read!(bin_io,vec)
                f(id,vec)
            end
        end
    end

    function load_milk_binaries(ids_path)
        data_dict = Dict{String,Vector{Float32}}()
        stream_binary_file(ids_path) do id,vec
            data_dict[id] = copy(vec)
        end
        groupsize_dict = Dict{String,Int}()
        for line in eachline(groupsize_path(ids_path))
            id,size = split(line,",")
            groupsize_dict[id] = parse(Int,size)
        end
        return data_dict,groupsize_dict
    end

    function write_milk_binaries(data_dict,groupsize_dict,ids_path)
        open(ids_path,"w") do id_io
            open(binary_path(ids_path),"w") do binary_io
                open(groupsize_path(ids_path),"w") do sizes_io
                    for (id,vec) in data_dict
                        println(id_io,id)
                        write(binary_io,vec)
                        println(sizes_io,id,",",groupsize_dict[id])
                    end
                end
            end
        end
    end

    function convert_input_csv_to_binaries(csv_path,ids_path)
        nan_count = 0
        open(ids_path,"w") do id_io
            open(binary_path(ids_path),"w") do binary_io
                open(groupsize_path(ids_path),"w") do sizes_io
                    for line in eachline(csv_path)
                        entry = split(line,",")
                        vec = [val == "" ? NaN32 : parse(Float32,val) for val in entry[2:end]]
                        if any(isnan,vec)
                            nan_count += 1
                            continue
                        end
                        println(id_io,entry[1])
                        write(binary_io,vec)
                        println(sizes_io,entry[1],",",1)
                    end
                end
            end
        end
        if nan_count > 0
            @warn "Dropped $(nan_count) rows with missing values from $(csv_path)"
        end
    end

    function partition_input_file(input_path,label,invariant_args)
        partition_dir = joinpath(invariant_args["output-dir"],"$(label).split")
        mkpath(partition_dir)

        function get_partitioned_input_path(p)
            partition_label = "partition_$(lpad(string(p),8,'0'))"
            return joinpath(partition_dir,"$(label).$(partition_label).ids")
        end

        p = 1
        data_dict = Dict{String,Vector{Float32}}()
        groupsize_dict = Dict{String,Int}()
        open(groupsize_path(input_path)) do sizes_io
            stream_binary_file(input_path) do id,vec
                size_id,size = split(readline(sizes_io),",")
                if size_id != id
                    error("Row mismatch between $(input_path) and its sizes file: $(id) vs $(size_id)")
                end
                data_dict[id] = copy(vec)
                groupsize_dict[id] = parse(Int,size)
                if length(data_dict) == invariant_args["partition-size"]
                    write_milk_binaries(data_dict,groupsize_dict,get_partitioned_input_path(p))
                    empty!(data_dict)
                    empty!(groupsize_dict)
                    p += 1
                end
            end
        end
        if !isempty(data_dict)
            write_milk_binaries(data_dict,groupsize_dict,get_partitioned_input_path(p))
        end
        return partition_dir
    end

    function batch_partitioned_files(partition_dir,label,invariant_args)

        function get_batch_directory(b)
            batch_label = "batch_$(lpad(string(b),8,'0'))"
            return joinpath(partition_dir,"$(label).$(batch_label).work")
        end

        pattern = "*.ids"
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
                    mv(binary_path(path),binary_path(updated_path))
                    mv(groupsize_path(path),groupsize_path(updated_path))
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
                mv(binary_path(path),binary_path(updated_path))
                mv(groupsize_path(path),groupsize_path(updated_path))
            end
        end
        return files,batches
    end

    function partition_and_batch_input_files(input_path,label,invariant_args)
        partition_dir = partition_input_file(input_path,label,invariant_args)
        input_paths,batch_dirs = batch_partitioned_files(partition_dir,label,invariant_args)
        return partition_dir,input_paths,batch_dirs
    end


    function open_file_write(f::Function, path::AbstractString; gzip::Bool=true)
        stream = gzip ? GzipCompressorStream(open(path,"w")) : open(path,"w")
        try
            return f(stream)  # Pass the stream to the function
        finally
            close(stream)  # Ensure the stream is closed properly
        end
    end

    function open_file_read(f::Function, path::AbstractString; gzip::Bool=true)
        stream = gzip ? GzipDecompressorStream(open(path,"r")) : open(path, "r")
        try
            return f(stream)  # Pass the stream to the function
        finally
            close(stream)  # Ensure the stream is closed properly
        end
    end

    function clean_directory(work_dir,partition_dir,exclusion_set)
        rm(partition_dir,recursive=true)
        for path in glob("*.representatives.ids",work_dir)
            if path in exclusion_set continue end
            rm(path)
            rm(binary_path(path))
            rm(groupsize_path(path))
        end
        for path in glob("*.input.ids",work_dir)
            if islink(path) && !isfile(path) # broken symlink: its representatives were deleted above
                rm(path)
                rm(binary_path(path))
                rm(groupsize_path(path))
            end
        end
        return
    end

    function final_cleanup(output_dir)
        for path in [glob("*.ids",output_dir); glob("*.bin",output_dir); glob("*.sizes",output_dir)]
            rm(path)
        end
        return
    end

    function write_values_as_txt(values,path)
        open(path,"w") do handle
            for value in values
                write(handle,"$value\n")
            end
        end
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

    function attempt_to_load_cache(path)
        cache_dict = nothing
        if !isnothing(path) && isfile(path)
            cache_dict,_ = load_milk_binaries(path)
        end
        return cache_dict
    end

    function attempt_to_load_previous_groups(path)
        previous_groups = nothing
        if !isnothing(path) && isfile(path)
            previous_groups_dict = load_groups_as_dictionary(path)
            previous_groups = previous_groups_dict["groups"]
        end
        return previous_groups
    end

    function prepare_milk_input(csv_path)
        milk_input_dir = joinpath(dirname(csv_path),"milk_input")
        file_label = replace(basename(csv_path),r"\.csv$" => "")

        ids_path = joinpath(milk_input_dir,"$(file_label).ids")
        source_path = joinpath(milk_input_dir,"$(file_label).source")
        source = "$(filesize(csv_path)),$(mtime(csv_path))"
        if isfile(source_path) && read(source_path,String) == source
            @info "Using existing MILK binaries: $(ids_path)"
            return ids_path
        end

        @info "Converting $(csv_path) to MILK binaries: $(ids_path)"
        mkpath(milk_input_dir)
        tmp_ids_path = joinpath(milk_input_dir,"$(file_label).tmp_$(getpid()).ids")
        convert_input_csv_to_binaries(csv_path,tmp_ids_path)
        mv(tmp_ids_path,ids_path,force=true)
        mv(binary_path(tmp_ids_path),binary_path(ids_path),force=true)
        mv(groupsize_path(tmp_ids_path),groupsize_path(ids_path),force=true)
        write(source_path,source)
        return ids_path
    end

end