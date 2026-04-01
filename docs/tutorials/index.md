# Tutorials

This section provides more detailed examples applying MILK to actual datasets and some of the downstream analyses that can be facilitated.

```{toctree}
:maxdepth: 2

tutorial.pbmc3k_dataset.network_visualization.md
tutorial.mnist_dataset.downsampling_visualization.md
tutorial.tfatlas_dataset.capturing_tf_induced_changes.md
```

## Use cases

### 1. Generating information-rich subsamples

As the scale of datasets continues to grow, downsampling the data into computationally tractable sizes is a processing step that is becoming increasingly necessary. However, this step can effectively result in 'throwing away' almost the all of the data.

> For example, many downstream biological inference algorithms cannot scale beyond 10<sup>5-6</sup> cells. With the emergence of massive single-cell resources containing upwards of 10<sup>7-8</sup> cells, downsampling can result in the loss of 90-99.9% of information.

Class imbalance, a phenomenon in which some classes are significantly more abundant than other classes, can also introduce a multitude of biases and is non-trivial to address at scale.

MILK trees can be utilized to generate representative subsamples with or without the use of additional information (e.g., metadata information). Termed label-based or label-free subsampling, respectively, this can be described as the task of selecting a set of internal nodes (according to a target sample size) with no overlap (between leaves corresponding to each clade, as defined by its respective internal node) while simultaneously maximizing coverage over the global tree structure.

Because each subsampled cell can refer to a clade in the global MILK tree, this represents a principled framework to retain information in a data-driven manner.

> Using single-cell transcriptomic profiles as an example, explicit information retention can refer to the aggregation of gene expression profiles to construct "metacells", whereas implicit information retention can be the augmentation of subsamples in downstream analyses via weighting cells with the associated information.

### 2. Global mapping of subpopulation dynamics

The pairwise relationships encoded in the global MILK hierarchy represents a powerful framework to explain changes between individual objects. While many approaches can identify such changes, they typically operate on slices at scale, making it difficult to ascertain how meaningful these changes are. In contrast, MILK trees can encode complex interactions all in a cohesive data structure, which allows a more clear understanding of *how* different they are.

This can be particularly powerful to capture transitions of between states in an interpretable framework.

> An example of this can be the characterizing cellular perturbations, as reflected in their resulting gene expression profiles (e.g., CRIPSR-mediated perturbation of transcription factors to study cell differentiation).

### 3. Multi-resolution analyses

The MILK tree structure is inherently conducive to multi-resultion analyses. In addition, there are efficient algorithms designed to work with hierarchical encodings of datasets. Together, this makes MILK suitable to carry out more comprehensive analyses focusing on how signal emerges in both local and global contexts.

> An example can entail characterizing the transcriptomic signatures of cell types across many species in a single tree.

### 4. Hierarchical clustering

### 5. Outlier detection

### 6. Unbiased nearest-neighbor detection (arbitrary *k*)
