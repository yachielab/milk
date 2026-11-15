module IncrementalFileHandling

    using Logging

    using ..FileHandling: open_file_read,open_milk_files,close_milk_files,mv_milk_files,
                          load_milk_binaries,load_groups_as_dictionary,group_meta_path,metadata_path

    export prepare_incremental_milk_input

    function incremental_ids_path(base_path,stage)
        return "$(base_path).stage_$(lpad(stage,8,'0')).ids"
    end

    function prepare_incremental_milk_input(csv_path,incremental_path,label)
        milk_input_dir = joinpath(dirname(csv_path),"$(label).input")
        file_label = replace(basename(csv_path),r"\.csv(\.gz)?$" => "")
        mkpath(milk_input_dir)

        base_path = joinpath(milk_input_dir,"$(file_label)")
        stages_path = joinpath(milk_input_dir,"$(file_label).stages.tsv")
        source_path = joinpath(milk_input_dir,"$(file_label).source")
        source = "$(filesize(csv_path)),$(mtime(csv_path)),$(filesize(incremental_path)),$(mtime(incremental_path))"

        n_columns = 0
        if isfile(stages_path)
            first_stage_ids_path = String(split(readline(stages_path),'\t')[2])
            n_columns = length(split(readline(group_meta_path(first_stage_ids_path)),'\t'))
        end
        if isfile(source_path) && read(source_path,String) == source && n_columns == 6
            @info "Using existing stage-stratified MILK binaries: $(stages_path)"
            return stages_path
        end

        @info "Converting input file to stage-stratified (incremental) MILK binaries: $(stages_path)"
        @info "\tInput: $(csv_path)"
        @info "\tStages: $(incremental_path)"
        tmp_base_path = joinpath(milk_input_dir,"$(file_label).tmp_$(getpid())")
        stage_dict = read_incremental_ids_as_dictionary(incremental_path)
        stage_counts = convert_input_csv_to_stage_binaries(csv_path,tmp_base_path,stage_dict)
        open(stages_path,"w") do io
            for stage in sort(collect(keys(stage_counts)))
                stage_ids = incremental_ids_path(base_path,stage)
                mv_milk_files(incremental_ids_path(tmp_base_path,stage),stage_ids;force=true)
                println(io,join([stage,stage_ids,stage_counts[stage]],"\t"))
            end
        end
        write(source_path,source)
        return stages_path
    end

    function read_incremental_ids_as_dictionary(path)
        """
        id -> stage (integer, 1 = earliest). Stage values in the table can be any real numbers
        (e.g., 8.5, 12, -1); they are sorted and mapped to 1, 2, 3, ...
        """
        value_dict = Dict{String,Float64}()
        label_dict = Dict{Float64,String}() # value -> as written in the table (first occurrence), for logging
        for line in eachline(path)
            id,value = split(line,"\t")
            parsed_value = tryparse(Float64,strip(value))
            if isnothing(parsed_value) || !isfinite(parsed_value)
                error("Stage of $(id) in $(path) is not a real number: $(value)")
            end
            value_dict[id] = parsed_value
            if !haskey(label_dict,parsed_value)
                label_dict[parsed_value] = String(strip(value))
            end
        end
        if isempty(value_dict)
            error("No ids in $(path)")
        end

        count_dict = Dict{Float64,Int}()
        for value in values(value_dict)
            count_dict[value] = get(count_dict,value,0) + 1
        end

        rank_dict = Dict{Float64,Int}()
        @info "Stage mapping ($(path)): value in table -> stage (number of ids)"
        for (rank,value) in enumerate(sort(collect(keys(count_dict))))
            rank_dict[value] = rank
            @info "\t$(label_dict[value]) -> $(rank) ($(count_dict[value]) ids)"
        end

        stage_dict = Dict{String,Int}()
        for (id,value) in value_dict
            stage_dict[id] = rank_dict[value]
        end
        return stage_dict
    end

    function convert_input_csv_to_stage_binaries(csv_path,base_path,stage_dict)
        """ one set of binaries per stage: <base_path>.stage_XXXXXXXX.ids/.bin/.group.tsv; returns stage => number of rows """
        files = Dict{Int,Any}()
        stage_counts = Dict{Int,Int}()
        n_rows = 0
        n_rows_without_stage = 0
        n_missing_values = 0
        open_file_read(csv_path,gzip=endswith(csv_path,".gz")) do instream
            for line in eachline(instream)
                n_rows += 1
                entry = split(line,",")
                if !haskey(stage_dict,entry[1])
                    n_rows_without_stage += 1
                    continue
                end
                vec = [val == "" ? NaN32 : parse(Float32,val) for val in entry[2:end]]
                if any(isnan,vec)
                    n_missing_values += 1
                    continue
                end
                stage = stage_dict[entry[1]]
                if !haskey(files,stage)
                    files[stage] = open_milk_files(incremental_ids_path(base_path,stage),"w")
                    stage_counts[stage] = 0
                end
                println(files[stage].ids,entry[1])
                write(files[stage].bin,vec)
                println(files[stage].group,join([entry[1],1,0.0,0.0,stage,0],"\t")) # id, size, spread, radius, stage, n_cached
                stage_counts[stage] += 1
            end
        end
        for milk_files in values(files)
            close_milk_files(milk_files)
        end

        n_rows_with_stage = n_rows - n_rows_without_stage
        if n_rows_with_stage == 0
            error("No ids shared between $(csv_path) and the stage table")
        end
        if n_rows_without_stage > 0
            @warn "Dropped $(n_rows_without_stage) rows of $(csv_path) whose id is not in the stage table"
        end
        if length(stage_dict) > n_rows_with_stage
            @warn "$(length(stage_dict)-n_rows_with_stage) ids in the stage table are not in $(csv_path)"
        end
        if n_missing_values > 0
            @warn "Dropped $(n_missing_values) rows with missing values from $(csv_path)"
        end
        return stage_counts
    end

end