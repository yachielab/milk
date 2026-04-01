<img src="https://raw.githubusercontent.com/yachielab/milk/main/imgs/milk_logo_repo.png" width=800>

[![Documentation Status](https://readthedocs.org/projects/milk-trees/badge/?version=latest)](https://milk-trees.readthedocs.io/en/latest/)

# Paper

Stay tuned...

# Documentation

For more detailed information, please refer to the [documentation](https://milk-trees.readthedocs.io/en/latest/).

# Overivew

To address the rapid growth of high-dimensional data in biology and other fields, we present MILK, a scalable approach for capturing highly granular relationships between objects (e.g., single cells). By recursively stratifying the population into groups of highly similar objects based on pairwise distance, MILK builds a hierarchical organization of objects that preserves both global structure and rare local subpopulations. This provides a framework for multi-resolution analysis of ultra-large datasets in a tractable and interpretable manner.

<img src="https://raw.githubusercontent.com/yachielab/milk/main/imgs/similarity_search_diagram_repo.png">

The main process of MILK involves grouping highly similar objects using a data-driven similarity threshold defined by a lower-tail percentile (e.g., 0.1%) of distances. In a single-pass through the dataset, each candidate object is compared to the representatives of existing groups. If the candidate has a distance below the threshold to any group, it is mapped to the group it is most similar to; otherwise, it is considered distinct, in which case it forms a new group and becomes its representative. The first candidate object forms a new group by default.

<img src="https://raw.githubusercontent.com/yachielab/milk/main/imgs/milk_recursive_framework_repo.png">

To handle datasets that exceed available memory, MILK employs a recursive, out-of-core strategy that involves partitioning the population into tractable subsets and processing them with a shared global similarity threshold. After identifying groups within each subpopulation, representatives are optimized to be the medoid of each group. All groups are then aggregated into a global list and, if tractable, an additional grouping process is performed on the representatives to resolve potential overlap across partitions.

By carrying forward only representatives, this pipeline can be applied recursively to progressively merge groups according to the percentile-based similarity threshold. If allowed to run to completion, a global hierarchical structure is produced over the complete population.

# Installation

## Option 1: Pre-compiled binary (Recommended)

Users can run MILK without installing Julia on your system.

> [!IMPORTANT]
> Pre-release note: currently, the x86 (64-bit) build is provided1

Download the latest `milk` build from the [Latest Releases](https://github.com/yachielab/milk/releases) page.

Extract binary:
```
tar -xzf milk-linux-x86_64.tar.gz
```

Export to `PATH` (replace `/path/to` with the actual download location):
```
export PATH="$PATH:/path/to/milk/bin"
```

## Option 2: MILK installation as a Julia module

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
## Verify successful installation
Confirm successful installation with the following call:
```
milk -h
```