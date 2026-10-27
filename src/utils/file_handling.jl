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
        metadata_path,
        binary_path,
        group_meta_path,
        stream_binary_file,
        load_milk_binaries,
        write_milk_binaries,
        convert_input_csv_to_binaries,
        prepare_milk_input,
        GroupCache,
        cache_meta_path,
        cache_binary_path,
        load_milk_cache,
        count_rows,
        load_group_spreads,
        mv_milk_files


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

    function group_meta_path(ids_path)
        if !endswith(ids_path,".ids")
            error("Expected a .ids path: $(ids_path)")
        end
        return replace(ids_path,r"\.ids$" => ".group.tsv")
    end

    function cache_meta_path(ids_path)
        return replace(ids_path,r"\.ids$" => ".cache.tsv")
    end

    function cache_binary_path(ids_path)
        return replace(ids_path,r"\.ids$" => ".cache.bin")
    end
    ####################################

    function count_rows(ids_path)
        if !isfile(cache_meta_path(ids_path))
            return countlines(ids_path)
        end
        return countlines(cache_meta_path(ids_path))
    end

    function load_groups_as_dictionary(path)

        threshold = split(readline(metadata_path(path)),'\t')[4]

        groups = Dict{String,Vector{String}}()
        spread_dict = Dict{String,Any}()
        specificity_dict = Dict{String,Any}()
        resolution_dict = Dict{String,Any}()
        radius_dict = Dict{String,Float64}()
        distances_dict = Dict{String,Vector{Float64}}()
        open_file_read(path,gzip=true) do file
            for line in eachline(file)
                info = JSON.parse(line)
                representative_id = info["representative_id"]
                groups[representative_id] = info["group"]
                spread_dict[representative_id] = isempty(info["distances"]) ? nothing : mean(info["distances"]) 
                specificity_dict[representative_id] = isempty(info["specificity"]) ? nothing : mean(info["specificity"]) 
                resolution_dict[representative_id] = length(info["distances"])
                radius_dict[representative_id] = isempty(info["distances"]) ? 0.0 : maximum(info["distances"])
                distances_dict[representative_id] = Float64.(info["distances"])
            end
        end
        groups_dict = Dict(
            "groups" => groups,
            "threshold" => threshold,
            "spread" => spread_dict,
            "specificity" => specificity_dict,
            "resolution" => resolution_dict,
            "radius" => radius_dict,
            "distances" => distances_dict
        )
        return groups_dict
    end

    function write_group_metadata(;path,label,stage,n_cached,threshold,n_input_objects,n_groups,n_comparisons)
        open(metadata_path(path),"w") do file
            # columns: label, stage, n_cached, threshold, n_input_objects, n_groups, n_comparisons
            fields = [label,stage,n_cached,threshold,n_input_objects,n_groups,n_comparisons]
            println(file,join(fields,'\t'))
        end
    end

    function write_group_results(;path,label,stage,n_cached,n_input_objects,n_groups,groups,direct_groupsize_dict,
                                  total_groupsize_dict,distances_dict,specificity_dict,threshold,n_comparisons,cache_dict=nothing)

        write_group_metadata(
            path=path,
            label=label,
            stage=stage,
            n_cached=n_cached,
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
                    "direct_group_size" => length(group),
                    "total_group_size" => total_groupsize_dict[representative_id],
                    "distances" => distances_dict[representative_id],
                    "specificity" => specificity_dict[representative_id]
                )
                if !isnothing(cache_dict)
                    cache = cache_dict[representative_id]
                    group_info["cache_ids"] = cache.cell_ids
                    group_info["cache_sizes"] = cache.subgroup_sizes
                    group_info["cache_radii"] = cache.subgroup_radii
                end
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
        for line in eachline(group_meta_path(ids_path))
            id,size = split(line,'\t')[1:2]
            groupsize_dict[id] = parse(Int,size)
        end
        return data_dict,groupsize_dict
    end

    function load_group_spreads(ids_path)
        """ group id => spread (mean squared distance of the group's cells to its centroid) """
        spread_dict = Dict{String,Float64}()
        for line in eachline(group_meta_path(ids_path))
            fields = split(line,'\t')
            spread_dict[String(fields[1])] = parse(Float64,fields[3])
        end
        return spread_dict
    end

    function write_milk_binaries(data_dict,groupsize_dict,ids_path;
                                 spread_dict=nothing,cache_dict=nothing)
        cache_flag = !isnothing(cache_dict)
        cache_meta_io = cache_flag ? open(cache_meta_path(ids_path),"w") : nothing
        cache_bin_io = cache_flag ? open(cache_binary_path(ids_path),"w") : nothing
        open(ids_path,"w") do id_io
            open(binary_path(ids_path),"w") do binary_io
                open(group_meta_path(ids_path),"w") do group_meta_io
                    for (id,vec) in data_dict
                        spread = isnothing(spread_dict) ? 0.0 : spread_dict[id]
                        println(id_io,id)
                        write(binary_io,vec)
                        if cache_flag
                            cache = cache_dict[id]
                            fields = [id,groupsize_dict[id],spread,length(cache.cell_ids)]
                            println(group_meta_io,join(fields,"\t"))
                            for k in eachindex(cache.cell_ids)
                                println(cache_meta_io,cache.cell_ids[k],'\t',cache.subgroup_sizes[k],'\t',cache.subgroup_radii[k])
                                write(cache_bin_io,cache.vecs[k])
                            end
                        else
                            fields = [id,groupsize_dict[id],spread]
                            println(group_meta_io,join(fields,"\t"))
                        end
                    end
                end
            end
        end
        if cache_flag
            close(cache_meta_io)
            close(cache_bin_io)
        end
    end

    function _read_group_cache(cache_meta_io,cache_bin_io,n_cached,d)
        """ reads the next n_cached cached cells (one group) from open cache files
        """
        cache = GroupCache(String[],Vector{Float32}[],Int[],Float64[])
        for _ in 1:n_cached
            cell_id,subgroup_size,subgroup_radius = split(readline(cache_meta_io),'\t')
            vec = Vector{Float32}(undef,d)
            read!(cache_bin_io,vec)
            push!(cache.cell_ids,String(cell_id))
            push!(cache.vecs,vec)
            push!(cache.subgroup_sizes,parse(Int,subgroup_size))
            push!(cache.subgroup_radii,parse(Float64,subgroup_radius))
        end
        return cache
    end

    function load_milk_cache(ids_path,data_dict,groupsize_dict;keep_ids=nothing)
        """ group id => GroupCache. Without cache files, each group in data_dict is cached as its own representative.
        With keep_ids, only those groups are loaded; the rest of the file is streamed past (memory bounded by keep_ids).
        """
        cache_dict = Dict{String,GroupCache}()
        if !isfile(cache_meta_path(ids_path))
            for (id,vec) in data_dict
                # representative caching
                cache_dict[id] = GroupCache([id],[vec],[groupsize_dict[id]],[0.0])
            end
            return cache_dict
        end
        d = isempty(data_dict) ? 0 : length(first(values(data_dict)))
        open(cache_meta_path(ids_path)) do cache_meta_io
            open(cache_binary_path(ids_path)) do cache_bin_io
                for line in eachline(group_meta_path(ids_path))
                    id,_,_,n_cached = split(line,'\t')
                    n_cached = parse(Int,n_cached)
                    if isnothing(keep_ids) || String(id) in keep_ids
                        cache_dict[String(id)] = _read_group_cache(cache_meta_io,cache_bin_io,n_cached,d)
                    else
                        for _ in 1:n_cached
                            readline(cache_meta_io)
                        end
                        skip(cache_bin_io,4*d*n_cached) # Float32
                    end
                end
            end
        end
        return cache_dict
    end

    function convert_input_csv_to_binaries(csv_path,ids_path)
        nan_count = 0
        open(ids_path,"w") do id_io
            open(binary_path(ids_path),"w") do binary_io
                open(group_meta_path(ids_path),"w") do group_meta_io
                    open_file_read(csv_path,gzip=endswith(csv_path,".gz")) do instream
                        for line in eachline(instream)
                            entry = split(line,",")
                            vec = [val == "" ? NaN32 : parse(Float32,val) for val in entry[2:end]]
                            if any(isnan,vec)
                                nan_count += 1
                                continue
                            end
                            println(id_io,entry[1])
                            write(binary_io,vec)
                            fields = [entry[1],1,0.0]
                            println(group_meta_io,join(fields,"\t"))
                        end
                    end
                end
            end
        end
        if nan_count > 0
            @warn "Dropped $(nan_count) rows with missing values from $(csv_path)"
        end
    end

    struct GroupCache
        cell_ids::Vector{String}
        vecs::Vector{Vector{Float32}}
        subgroup_sizes::Vector{Int}
        subgroup_radii::Vector{Float64}
    end

    function partition_input_file(input_path,label,invariant_args)
        partition_dir = joinpath(invariant_args["output-dir"],"$(label).split")
        mkpath(partition_dir)

        function _get_partitioned_input_path(p)
            partition_label = "partition_$(lpad(string(p),8,'0'))"
            return joinpath(partition_dir,"$(label).$(partition_label).ids")
        end

        function _write_partition()
            path = _get_partitioned_input_path(p)

            write_milk_binaries(data_dict,groupsize_dict,path;spread_dict=spread_dict,
                                cache_dict=(cache_flag ? cache_dict : nothing))
            empty!(data_dict)
            empty!(groupsize_dict)
            empty!(spread_dict)
            empty!(cache_dict)
            n_rows = 0
            p += 1
        end

        cache_flag = isfile(cache_meta_path(input_path))
        cache_meta_io = cache_flag ? open(cache_meta_path(input_path)) : nothing
        cache_bin_io = cache_flag ? open(cache_binary_path(input_path)) : nothing

        p = 1
        data_dict = Dict{String,Vector{Float32}}()
        groupsize_dict = Dict{String,Int}()
        spread_dict = Dict{String,Float64}()
        cache_dict = Dict{String,GroupCache}()
        n_rows = 0

        open(group_meta_path(input_path)) do group_meta_io
            stream_binary_file(input_path) do id,vec
                fields = split(readline(group_meta_io),'\t')
                if fields[1] != id
                    error("Row mismatch between $(input_path) and its group file: $(id) vs $(fields[1])")
                end
                data_dict[id] = copy(vec)
                groupsize_dict[id] = parse(Int,fields[2])
                spread_dict[id] = parse(Float64,fields[3])

                if cache_flag
                    n_cached = parse(Int,fields[4])
                    cache_dict[id] = _read_group_cache(cache_meta_io,cache_bin_io,n_cached,length(vec))
                    n_rows += n_cached
                else
                    n_rows += 1
                end

                if n_rows >= invariant_args["partition-size"]
                    _write_partition()
                end
            end
        end
        if !isempty(data_dict)
            _write_partition()
        end
        if cache_flag
            close(cache_meta_io)
            close(cache_bin_io)
        end
        return partition_dir
    end

    function mv_milk_files(src_ids_path,dst_ids_path;force=false)
        """ moves .ids/.bin/.group.tsv, plus .cache.tsv/.cache.bin if present """
        mv(src_ids_path,dst_ids_path;force=force)
        mv(binary_path(src_ids_path),binary_path(dst_ids_path);force=force)
        mv(group_meta_path(src_ids_path),group_meta_path(dst_ids_path);force=force)
        if isfile(cache_meta_path(src_ids_path))
            mv(cache_meta_path(src_ids_path),cache_meta_path(dst_ids_path);force=force)
            mv(cache_binary_path(src_ids_path),cache_binary_path(dst_ids_path);force=force)
        end
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
                    mv_milk_files(path,updated_path)
                end
            end
        else
            batch_dir = get_batch_directory(0)
            mkpath(batch_dir)
            push!(batches,batch_dir)
            for path in paths
                updated_path = joinpath(batch_dir,basename(path))
                push!(files,updated_path)
                mv_milk_files(path,updated_path)
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
            rm(group_meta_path(path))
            rm(cache_meta_path(path),force=true)
            rm(cache_binary_path(path),force=true)
        end
        for path in glob("*.input.ids",work_dir)
            if islink(path) && !isfile(path) # broken symlink: its representatives were deleted above
                rm(path)
                rm(binary_path(path))
                rm(group_meta_path(path))
                rm(cache_meta_path(path),force=true)
                rm(cache_binary_path(path),force=true)
            end
        end
        return
    end

    function final_cleanup(output_dir)
        for path in [glob("*.ids",output_dir); glob("*.bin",output_dir); glob("*.group.tsv",output_dir); glob("*.cache.tsv",output_dir)]
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

    function prepare_milk_input(csv_path,label)
        milk_input_dir = joinpath(dirname(csv_path),"$(label).input")
        file_label = replace(basename(csv_path),r"\.csv(\.gz)?$" => "")

        ids_path = joinpath(milk_input_dir,"$(file_label).ids")
        source_path = joinpath(milk_input_dir,"$(file_label).source")
        source = "$(filesize(csv_path)),$(mtime(csv_path))"

        if isfile(source_path) && read(source_path,String) == source && isfile(group_meta_path(ids_path)) &&
           length(split(readline(group_meta_path(ids_path)),'\t')) == 3
            @info "Using existing MILK binaries: $(ids_path)"
            return ids_path
        end

        @info "Converting input file to MILK binaries: $(ids_path)"
        @info "\tInput: $(csv_path)"
        @info "\tBinaries: $(milk_input_dir)"
        mkpath(milk_input_dir)
        tmp_ids_path = joinpath(milk_input_dir,"$(file_label).tmp_$(getpid()).ids")
        convert_input_csv_to_binaries(csv_path,tmp_ids_path)
        mv_milk_files(tmp_ids_path,ids_path;force=true)
        write(source_path,source)
        return ids_path
    end

end