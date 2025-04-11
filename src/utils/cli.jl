using Logging

function instantiate_invariant_args(args)
    invariant_args = Dict(
        "batch-size"       => args["batch-size"],
        "sample-size"      => args["sample-size"],
        "cache-size-limit" => args["cache-size-limit"],
        "partition-size"   => args["partition-size"],
        "merge-threshold"  => args["merge-threshold"],
        "percentile"       => args["percentile"],
        "metric"           => args["metric"],
        "seed"             => args["seed"],
        "threads"          => args["threads"],
        "job-scheduler"    => args["job-scheduler"],
        "job-account"      => args["job-account"],
        "job-name"         => args["job-name"],
        "job-memory"       => args["job-memory"],
        "environment-path" => args["environment-path"],
        "output-dir"       => abspath(args["output-dir"])
    )
    return invariant_args
end

function log_args(args)
    @info "Parsed Arguments:"
    for (key,value) in args
        @info "$(key): $(value)"
    end
    flush(stdout)
end


