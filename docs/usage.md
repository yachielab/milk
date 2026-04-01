# Usage

```
usage: milk -i INPUT-PATH [-n SAMPLE-SIZE] [-c CACHE-SIZE-LIMIT]
            [-p PARTITION-SIZE] [-M MERGE-THRESHOLD] [-m METRIC]
            [-l LABEL] [-t PERCENTILE] [-T THREADS] [-r SEED]
            [--verbose] [--force-overwrite] [--skip-reconstruction]
            [-o OUTPUT-DIR] [--hpc-mode] [-b BATCH-SIZE]
            [--job-scheduler JOB-SCHEDULER]
            [--job-account JOB-ACCOUNT] [--job-name JOB-NAME]
            [--job-memory JOB-MEMORY] [--job-time JOB-TIME]
            [--environment-path ENVIRONMENT-PATH]

A command-line tool to capture hierarchical relationships at scale.

optional arguments:
  -h, --help            show this help message and exit

Main arguments:
  -i, --input-path INPUT-PATH
                        Uncompressed CSV file. First column
                        corresponds to object IDs (e.g., cell
                        barcode). No header
  -n, --sample-size SAMPLE-SIZE
                        Sample size threshold for the number of groups
                        required to stop the recursive grouping
                        process. Recursion stops when the number of
                        (representative) objects is less than or equal
                        to this value. (type: Int64, default: 1)
  -c, --cache-size-limit CACHE-SIZE-LIMIT
                        How many (representative) objects to cache for
                        medoid optimization in subsequent recursions.
                        (type: Int64, default: 50000)
  -p, --partition-size PARTITION-SIZE
                        How many objects per partitioned file. (type:
                        Int64, default: 10000)
  -M, --merge-threshold MERGE-THRESHOLD
                        Threshold of when to carry out an additional
                        merging stratification process (if exceeded,
                        basic concatenation to join partitioned
                        groups) (type: Int64, default: 100000)
  -m, --metric METRIC   Distance metric for pairwise comparisons.
                        Permissible values:
                        {'cosine','euclidean','manhattan','hamming','jaccard','correlation'}
                        (default: "euclidean")
  -l, --label LABEL     Specify string label for file names.
  -t, --percentile PERCENTILE
                        Threshold for stratification process to group
                        cells. (type: Float64, default: 1.0)
  -T, --threads THREADS
                        Number of distributed processes for computing
                        pairwise comparisons. (TODO change to n_cpus)
                        (type: Int64, default: 1)
  -r, --seed SEED       Random seed. (type: Int64, default: 21)
  --verbose             Amount of information written to standard
                        out/err.
  --force-overwrite     Overwrite output directory if it already
                        exists.
  --skip-reconstruction
                        Only perform hierarchical
                        grouping/downsampling phase and NOT the
                        additional step of reconstructing vertices and
                        edges table for graph structure.
  -o, --output-dir OUTPUT-DIR
                        Directory to write intermediate and output
                        files. (default: "./milk.out")

[High-Performance Computing mode] Cluster computing for large-scale MILK runs:
  --hpc-mode            Activate HPC mode for MILK execution. If this
                        flag is not specified, all [HPC mode]
                        arguments will be ignored.
  -b, --batch-size BATCH-SIZE
                        [HPC mode] The number of partitioned input
                        files to be assigned per job. (type: Int64,
                        default: 50)
  --job-scheduler JOB-SCHEDULER
                        [HPC mode] Job scheduler to distribute
                        partitioned batches of stratification
                        processes. Supported schedulers:
                        {'slurm','sge'}
  --job-account JOB-ACCOUNT
                        [HPC mode] Account name associated with
                        account (only for SLURM job scheduler).
  --job-name JOB-NAME   [HPC mode] Job names (default: "milk")
  --job-memory JOB-MEMORY
                        [HPC mode] Memory of jobs (in gigabytes).
                        (type: Int64, default: 4)
  --job-time JOB-TIME   [HPC mode] Allotted time for jobs (only for
                        SLURM job scheduler). (default: "24:00:00")
  --environment-path ENVIRONMENT-PATH
                        [HPC mode] Bash script to set up environment
                        in submitted jobs.
```

## Example usage

Below is the minimal `milk` call in terminal:
```
milk -i input.csv
```

## Input

The input file is an uncompressed CSV containing no headers. The first column denotes object IDs, which will be interpreted as strings. Subsequent columns contain the aligned high-dimensional values associated with each object (row).

> Permissible values can refer to real numbers including gene expression counts, lower-dimensional embeddings (e.g., principal components, latent embeddings, etc...).

Below is an excerpt from the PBMC 3k dataset used in tutorial that was obtained with `head -n 5 ./tutorial/input.csv`:
```
AAACATACAACCAC-1,5.556233,0.2577139,-0.18681023,...
AAACATTGAGCTAC-1,7.20953,7.4819846,0.16270587,...
AAACATTGATCAGC-1,2.6944375,-1.5836583,-0.6631259,...
AAACCGTGCTTCCG-1,-10.143295,-1.3685299,1.2098124,...
AAACCGTGTATGCG-1,-1.1128161,-8.152788,1.3324049,...
...
```
> Currently, objects containing any missing data will be filtered.

## Output

MILK outputs a hierarchical representation of the input dataset as a tabular edge and node list format.

`edges.csv.gz` contains links from parent to child node/vertex. The source/target directionality is oriented so that the root (final iteration) of the tree is the source and the leaves (initial population) are the targets.
```
source,target
I1,CAATCGGAGAAACA-1
I1,TTCTACGAACGTAC-1
I2,CCTCGAACACTTTC-1
I3,GAACCTGATGAACC-1
I4,GAGGTTTGTAAGCC-1
```

`vertices.csv.gz` contains the following information for every (row) node/vertex:

* `node_id`: Unique ID for each node in the MILK tree (leaf IDs correspond to object IDs, whereas internal nodes take on an `I*` ID).
* `representative_id`: The object ID of the representative of the group. For internal nodes, this will correspond to the medoid object of that group.
* `group_size`: Clade size of the node/vertex. For any given internal node, its leaves correspond to the group members.
* `iteration`: Recursive iteration at which each grouping occurred.
* `threshold`: The percentile-based threshold at which each grouping occurred.
* `spread`: Approximate mean distance of each group member to its medoid.
* `specificity`: Approximate mean number of groups each member could viably be mapped to given the threshold at that iteration.
* `resolution`: The number of comparisons used to determine `spread` and `specificty`. This is based on the cache size. Greater the resolution, the more robust the approximation is likely to be.

## MILK heuristics

To achieve scalability, MILK implements heuristics to circumvent exhaustive calculation of pairwise comparisons, which scales quadratically with the number of objects.

Firstly, in the grouping process, each candidate object is only compared to the representatives of currently existing groups. The similarity threshold on the lower-extreme tail of the pairwise distribution ensures that only highly similar objects to the group representatives will merge.

Secondly, the partitioning of data into computationally tractable subsets dictates which objects can be potentially grouped. As a consequence, globally optimal groupings may not be identified. This is particularly relevant when the total number of groups across all partitions is intractable (i.e., above the merge threshold). To mitigate this MILK applies a shared global similarity threshold across all partitions. Related to this, the order of objects is preserved throughout the entire MILK execution process. This preservation can be leveraged to impart "prior knowledge" based on metadata information (e.g., pre-sort by cell type labels) or allow a more structure shuffle randomization of objects across multiple trials of MILK.

## High Performance Computing mode

For very large-scale datasets (i.e., on the order of tens to hundreds of millions of objects), MILK can be executed on HPC clusters. Given that partitioning of the entire population into tractable subsets becomes necessary with MILK at increasing dataset scales, parallelization by distributed computing processes (`-T` argument) becomes insufficient once there are thousands of partitions to process.

> For example, if you set a partition size of 10k cells (~50M pairwise comparisons), a dataset with 10M cells would generate 1k partitions.

The strategy to overcome this computational demand is to create batches of partitions (e.g., 50 partitions per batch), where each batch is submitted as as independent job on a HPC cluster. 

> The computing costs of 1k partitions can then be handled across 1000/50=20 jobs.

Batch size can be set by specifying the `-b` argument.

Note that if the `--hpc-mode` flag is not specified in the `milk` call, then any HPC-specific arguments will be ignored.

An example call is provided below:
```
milk \
-i input.csv \
-T 4 \
--hpc-mode \
--environment-path env_setup.sh \
--job-scheduler 'sge' \
--job-memory 8
```
The `--environment-path` argument expects the file path to a shell script that carrious out the activation of environment variables (e.g., conda environments), or exporting the milk executable to `PATH`. It will be invoked at the start of submitted job scripts as follows: `source env_setup.sh`.