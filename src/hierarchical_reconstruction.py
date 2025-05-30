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
    # parser.add_argument('-i','--input-path',type=str,required=True,help="")
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

def extract_iteration_label(path):
    iteration = os.path.basename(path).replace(".gz","").replace(".groups.jsonl","")
    iteration = int(iteration.split(".")[-1].split("_")[-1])
    return iteration

def open_file_read(path):
    if path.endswith(".gz"):
        return gzip.open(path,"rt")
    else:
        return open(path,"r")

def hierarchical_reconstruction():
    args = parse_arguments()

    clades_dict = {}
    ancestor_map = {}
    # with gzip.open(args.input_path,"rt") as handle:
    #     for line in handle:
    #         representative_id = line.strip().split(",")[0]
    #         clades_dict[representative_id] = instantiate_leaf_node(representative_id)
    # n = len(clades_dict)

    # pathlist = sorted([
    #     p for p in glob(os.path.join(args.out_dir,"*.groups.jsonl.gz"))
    #     if not (p.endswith(".concatenated.groups.jsonl.gz") or p.endswith(".merged.groups.jsonl.gz"))
    # ])

    pathlist = []
    pattern = os.path.join(args.out_dir,"*.groups.jsonl*")
    for path in sorted(glob(pattern)):
        field = os.path.basename(path).replace(".gz","").replace(".groups.jsonl","").split(".")[-1]
        if (field == "merged" or field == "concatenated"): continue
        pathlist.append(path)


    print("Instantiating clades_dict on initial recursion (iteration 1)",flush=True)
    groups_path = pathlist[0]
    iteration = extract_iteration_label(groups_path)
    n = 0
    # with gzip.open(groups_path,"rt") as handle:
    with open_file_read(groups_path) as handle:
        for line in handle:
            group_dict = json.loads(line)
            clade = instantiate_node_from_group_dictionary(group_dict,iteration)
            for sample_id in group_dict["group"]:
                clade.clades.append(instantiate_leaf_node(sample_id))
                ancestor_map[sample_id] = clade.name
                n += 1
            clades_dict[clade.name] = clade

    print(f"Reconstructing hierarchical groupings ({len(pathlist)} iterations)")
    for groups_path in pathlist[1:]:
        iteration = extract_iteration_label(groups_path)
        print(f"\tIteration: {iteration+1}",flush=True)
        # with gzip.open(groups_path,"rt") as handle:
        with open_file_read(groups_path) as handle:
            for line in handle:
                group_dict = json.loads(line)
                ancestors = {ancestor_map[sample_id] for sample_id in group_dict["group"]}
                clade = instantiate_node_from_group_dictionary(group_dict,iteration)
                clade.clades = [clades_dict.pop(sample_id) for sample_id in ancestors]
                clades_dict[clade.name] = clade

                for sample_id in group_dict["group"]:
                    ancestor_map[sample_id] = clade.name
                
        assert np.sum([c.count_terminals() for c in clades_dict.values()]) == n
    
    assert len(clades_dict) == 1
    tree = Tree(clades_dict[list(clades_dict)[0]])

    print("Finished hierarchical reconstruction!",flush=True)
    print("Creating vertices and edges mappings...",flush=True)

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
            id_map[label] = f"I{len(id_map)}"

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
                id_map[sublabel] = f"I{len(id_map)}"

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
