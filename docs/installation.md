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