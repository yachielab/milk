<img src="https://raw.githubusercontent.com/yachielab/milk/main/imgs/milk_logo_repo.png" width=800>

# Paper

TODO

# Overivew

To address the rapid growth of high-dimensional data in biology and other fields, we present MILK, a scalable approach for capturing highly granular relationships between objects (e.g., single cells). By recursively stratifying the population into groups of highly similar objects based on pairwise distance, MILK builds a hierarchical organization of objects that preserves both global structure and rare local subpopulations. This provides a framework for multi-resolution analysis of ultra-large datasets in a tractable and interpretable manner.

<img src="https://raw.githubusercontent.com/yachielab/milk/main/imgs/similarity_search_diagram_repo.png">

The main process of MILK involves grouping highly similar objects using a data-driven similarity threshold defined by a lower-tail percentile (e.g., 0.1%) of distances. In a single-pass through the dataset, each candidate object is compared to the representatives of existing groups. If the candidate has a distance below the threshold to any group, it is mapped to the group it is most similar to; otherwise, it is considered distinct, in which case it forms a new group and becomes its representative. The first candidate object forms a new group by default.

<img src="https://raw.githubusercontent.com/yachielab/milk/main/imgs/milk_recursive_framework_repo.png">

To handle datasets that exceed available memory, MILK employs a recursive, out-of-core strategy that involves partitioning the population into tractable subsets and processing them with a shared global similarity threshold. After identifying groups within each subpopulation, representatives are optimized to be the medoid of each group. All groups are then aggregated into a global list and, if tractable, an additional grouping process is performed on the representatives to resolve potential overlap across partitions.

By carrying forward only representatives, this pipeline can be applied recursively to progressively merge groups according to the percentile-based similarity threshold. If allowed to run to completion, a global hierarchical structure is produced over the complete population.


# Installation

Install Julia from the official [website](https://julialang.org/downloads/).

> Note that MILK was developed and tested using Julia version 1.11.3.
---
Install MILK using Julia's package manager:

In REPL:
```
using Pkg
Pkg.add(url="https://github.com/yachielab/milk.git")
```

Or, from the shell:
```
julia -e 'using Pkg; Pkg.add(url="https://github.com/yachielab/milk.git")'
```
---
Optionally, users can export the source directory to `PATH`:
```
echo 'export PATH="$PATH:./milk"
```
---
Confirm successful installation with the following call:
```
milk -h
```

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

Below is the minimal `milk` call :
```
milk -i input.csv
```

### Input

The input file is an uncompressed CSV containing no headers. The first column denotes object IDs, which will be interpreted as strings. Subsequent columns contain the aligned high-dimensional values associated with each object.

Permissible values can refer to real numbers including gene expression counts, lower-dimensional embeddings (e.g., principal components, latent embeddings, etc...).

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

### Output

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
