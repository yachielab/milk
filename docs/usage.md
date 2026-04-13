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

Capturing the global landscape of populations typically requires computation of an all-all pairwise matrix, where the number of comparisons scales quadratically with the number of individuals. This becomes intractable when datasets contain even tens to hundreds of thousands of individuals.

MILK implements a couple of heuristics to circumvent exhaustive calculation of pairwise comparisons.

First, in the grouping process, each candidate object is only compared to the representatives of currently existing groups. While this significantly reduces the number of comparisons, it heavily relies on the representatives being *good*. However, by applying a similarity threshold on the lower-extreme tail of the pairwise distribution, MILK aims to maximize resolution, ensuring only highly similar objects to the group representatives will merge (e.g., exhibiting similarity to the 99.9th percentile).

Partitioning of data into computationally tractable subsets is another heuristic that dictates which objects can be grouped. At scale, the likely consequence is that globally optimal groupings will not be identified. This is particularly relevant when the total number of groups across all partitions is intractable (i.e., above the merge threshold). MILK applies a shared global similarity threshold across all partitions, which should practically identify an upper bound to the number of possible groups at a given iteration. Related to this, the order of objects is preserved throughout the entire MILK execution process. This preservation can be leveraged to impart "prior knowledge" based on metadata information (e.g., pre-sort by cell type labels) or allow a more structure shuffle randomization of objects across multiple trials of MILK (planned extension).

Calculation of pairwise similarity/distance also depends on the dimensionality of the data. In the context of high-dimensional biological measurements (e.g., transcriptomic profiles), we typically recommend applying MILK to lower-dimensional embeddings of cells (e.g., principal components analysis, multi-dimensional scaling, non-negative matrix factorization, or latent embeddings from machine learning models). In addition to significantly speeding up computing speeds, it escapes from having to explicitly handle technical noise and sparsity (e.g., dropout) in pairwise comparisons -- this is something that should to be addressed but is out of the scope of MILK.

Overall, while suboptimal, these decisions have shown practical ability to capture biologically meaningful relationships at scale. Given the fact that single-cell measurements are riddled with technical (dropout, batch effect, sampling bias), biological (stochasticity, transcriptional bursting), and computational biases (processing, methodological assumptions), exact preservation of pairwise relationships is often neither achievable nor necessary at the scale of millions to hundreds of millions of cells. In this setting, approximate representations that preserve global structure can remain highly informative (e.g., identify emergent properties), particularly when integrating across large and heterogeneous datasets.

## High Performance Computing mode

The efficiency of the MILK algorithm is contingent on parallelization of computation on the tractable partitions of data. This is predominantly based on the number of CPUs that can be utilized during a MILK run (`-T` argument).

For very large-scale datasets (i.e., on the order of tens to hundreds of millions of objects), MILK can be executed on HPC clusters. Given that partitioning of the entire population into tractable subsets becomes necessary with MILK at increasing dataset scales, parallelization by 
distributed computing processes becomes insufficient once there are thousands of partitions to process.

> For example, if you set a partition size of 10k cells (~50M pairwise comparisons), a dataset with 10M cells would generate 1k partitions.

The strategy to overcome this computational demand is to create batches of partitions (e.g., 50 partitions per batch), where each batch is submitted as as independent job on a HPC cluster. 

> The computing costs of 1k partitions can then be handled across 1000/50=20 jobs.

Batch size can be set by specifying the `-b` argument. Users should set this argument with the job allocation time in mind -- if job allocation is slower then it should be much higher to prevent it from becoming the bottleneck.

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