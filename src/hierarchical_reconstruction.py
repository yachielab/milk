import os
import gzip
import json

import numpy as np
import pandas as pd
from glob import glob
from pathlib import Path
from collections import defaultdict
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree,Clade

from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description="")
    parser.add_argument('-i','--input-path',type=str,required=True,help="")
    parser.add_argument('-l','--label',type=str,required=True,help="")
    parser.add_argument('-o','--out-dir',type=str,required=True,help="")
    return parser.parse_args()

def instantiate_leaf_node(representative_id):
    clade = Clade(name=representative_id)
    clade.threshold = 0
    clade.iteration = -1
    clade.distances = None
    clade.specificity = None
    clade.resolution = None
    return clade

def instantiate_node_from_group_dictionary(group_dict,iteration):
    clade = Clade(name=group_dict["representative_id"])
    clade.threshold = group_dict["threshold"]
    clade.iteration = iteration
    clade.distances = group_dict["distances"]
    clade.specificity = group_dict["specificity"]
    clade.resolution = None if (clade.distances is None) else len(group_dict["distances"])
    return clade

def hierarchical_reconstruction():
    args = parse_arguments()

    clades_dict = {}
    with gzip.open(args.input_path,"rt") as handle:
        for line in handle:
            representative_id = line.strip().split(",")[0]
            clades_dict[representative_id] = instantiate_leaf_node(representative_id)
    n = len(clades_dict)

    pathlist = sorted([
        p for p in glob(os.path.join(args.out_dir,"*.groups.jsonl.gz"))
        if not (p.endswith(".concatenated.groups.jsonl.gz") or p.endswith(".merged.groups.jsonl.gz"))
    ])
    print(f"Reconstructing hierarchical groupings ({len(pathlist)} iterations)")
    for groups_path in pathlist:
        i = int(os.path.basename(groups_path).replace(".groups.jsonl.gz","").split("_")[-1])
        # print(f"\tIteration: {i+1}",flush=True)
        with gzip.open(groups_path,"rt") as handle:
            for line in handle:
                group_dict = json.loads(line)
                clade = instantiate_node_from_group_dictionary(group_dict,i)
                clade.clades = [clades_dict.pop(x) for x in group_dict["group"] if x in clades_dict]
                clades_dict[clade.name] = clade

        assert np.sum([c.count_terminals() for c in clades_dict.values()]) == n
    
    assert len(clades_dict) == 1
    tree = Tree(clades_dict[list(clades_dict)[0]])

    tree_dict = defaultdict(list)
    id_map = {}
    metadata_dict = {
        "node_id":[],
        "representative_id":[],
        "group_size":[],
        "iteration":[],
        "threshold":[],
        "spread":[],
        "specificity":[],
        "resolution":[]
    }

    stack = [tree.clade]
    while stack:
        clade = stack.pop(0)

        label = f"{clade.name}|{clade.iteration+1}"
        if label not in id_map:
            id_map[label] = len(id_map)

        metadata_dict["node_id"].append(id_map[label])
        metadata_dict["representative_id"].append(clade.name)
        metadata_dict["group_size"].append(clade.count_terminals())
        metadata_dict["iteration"].append(clade.iteration+1)
        metadata_dict["threshold"].append(clade.threshold)

        spread      = np.mean(clade.distances) if clade.distances else None
        specificity = np.mean(clade.specificity) if clade.specificity else None
        resolution  = 0 if spread is None else len(clade.distances)
        metadata_dict["spread"].append(spread)
        metadata_dict["specificity"].append(specificity)
        metadata_dict["resolution"].append(resolution)

        clade.name = str(id_map[label])
        for subclade in clade:

            sublabel = f"{subclade.name}|{subclade.iteration+1}"
            if sublabel not in id_map:
                id_map[sublabel] = len(id_map)

            tree_dict[clade.name].append(id_map[sublabel])
            stack.append(subclade)

    metadata_df = pd.DataFrame(metadata_dict)

    tree_dir = os.path.join(args.out_dir,"tree_representation")
    Path(tree_dir).mkdir(exist_ok=True,parents=False)

    print("\tWriting output files...")

    newick_path = os.path.join(tree_dir,f"{args.label}.tree.nwk")
    print(f"\t\tNewick path: {newick_path}",flush=True)
    Phylo.write(tree,newick_path,"newick")

    vertices_path = os.path.join(tree_dir,f"{args.label}.vertices.csv.gz")
    print(f"\t\tVertices path: {vertices_path}",flush=True)
    metadata_df.to_csv(vertices_path,compression="gzip",index=False)

    edges_path = os.path.join(tree_dir,f"{args.label}.edges.csv.gz")
    print(f"\t\tEdges path: {edges_path}",flush=True)
    with gzip.open(edges_path,"wt") as handle:
        for node,group in tree_dict.items():
            for child_node in group:
                handle.write(f"{node},{child_node}\n")

if __name__ == "__main__":
    hierarchical_reconstruction()
