# CRPEC
## Cluster-Refine-Predict-Ensemble-Cluster
An unsupervised clustering algorithm for scRNA-seq

## Overview

### Initialization
#### Scripts Involved:
- initialize_directories.sh
- trident_preprocess_to_h5ad_barcodes.py
### Cluster
Initial clustering on smaller subsets to produce 'pseudolabels'
#### Scripts Involved:
- h5ad_to_h5seurat.R
- preprocess_sc3.R
### Refine
Refining of previously obtained pseudolabels via recursive permutation test approach.
#### Scripts Involved:
- refine_clusters.py & refine_clusters.sh
### Predict
Train supervised predictive model with obtained pseudolabels to infer clusters of whole dataset.
#### Scripts Involved:
- predict_clusters.py & predict_clusters.sh
### Ensemble
Generate consensus on predicted clusters of multiple runs 
#### Scripts Involved:
- get_consensus_final_cluster_dict.py
- get_ensemble_adatas.py
#### Visualization Tools:
- get_cluster_map.ipynb
### Cluster
Cluster inter-sample gene profiles to generate dendrogram 
- get_gene_profiles_dtf.py
#### Visualization Tools:

