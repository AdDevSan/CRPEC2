# CRPEC
## Cluster-Refine-Predict-Ensemble-Cluster
An unsupervised clustering algorithm for scRNA

## Overview
### Cluster
Initial clustering on smaller subsets to produce 'pseudolabels'
#### Scripts Involved:
- initialize_directories.sh
- trident_preprocess_to_h5ad_barcodes.py
- h5ad_to_h5seurat.R
- preprocess_sc3.R
### Refine
Refining of previously obtained pseudolabels
#### Scripts Involved:
- refine_clusters.sh
- refine_clusters.py
### Predict
Train supervised predictive model with pseudolabels of dataset subset to infer clusters of whole dataset.
#### Scripts Involved:
- predict_clusters.py

### Ensemble
Generate consensus on predicted clusters of multiple runs 

### Cluster
Cluster inter-sample gene profiles to generate dendrogram 

