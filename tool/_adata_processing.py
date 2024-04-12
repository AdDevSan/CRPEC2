import os
import glob
import pandas as pd
import scanpy as sc
from scipy.io import mmread
import numpy as np
import pickle
import json
def find_file(directory, file_suffix):
    """Find a file in a given directory with the specified suffix."""
    files = glob.glob(os.path.join(directory, f'*{file_suffix}'))
    if len(files) == 0:
        raise FileNotFoundError(f"No file with suffix '{file_suffix}' found in directory '{directory}'")
    return files[0]

def read_mtx_to_adata(matrix_file, genes_file, barcodes_file):
    """Read matrix, genes, and barcodes files into an AnnData object."""
    matrix = mmread(matrix_file).T.tocsr()  # Transpose to have genes in columns
    genes = pd.read_csv(genes_file, header=None, sep='\t')[1].tolist()  # Assuming second column contains gene names
    barcodes = pd.read_csv(barcodes_file, header=None, sep='\t')[0].tolist()
    adata = sc.AnnData(X=matrix, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
    return adata
def get_subdirectories(base_directory):
    subdirectories = [os.path.join(base_directory, d) for d in os.listdir(base_directory) if os.path.isdir(os.path.join(base_directory, d))]
    return subdirectories


# takes base directory containing subdirectories (containing: matrix, barcodes, features) and makes adata for each subdirectory, returns adatas list 
def get_adatas_from_base_directory(base_directory):
    subdirectories = get_subdirectories(base_directory)
    adatas = []
    for directory in subdirectories:
        # Find each file in the directory
        matrix_file = find_file(directory, 'matrix.mtx.gz')
        barcodes_file = find_file(directory, 'barcodes.tsv.gz')
        genes_file = find_file(directory, 'features.tsv.gz')
        
        # Load the data into an AnnData object
        adata = read_mtx_to_adata(matrix_file, genes_file, barcodes_file)
        
        # Ensure the uniqueness of variable and observation names within each AnnData object
        adata.var_names_make_unique()
        adata.obs_names_make_unique('-')  # Using '-' to append a unique suffix if needed
        
        # Append to the list of AnnData objects
        adatas.append(adata)

    return adatas

def get_adata_from_trident(trident_directory_path):
    directory = trident_directory_path
    matrix_file = find_file(directory, 'matrix.mtx.gz')
    barcodes_file = find_file(directory, 'barcodes.tsv.gz')
    genes_file = find_file(directory, 'features.tsv.gz')

    adata = read_mtx_to_adata(matrix_file, genes_file, barcodes_file)
    adata.var_names_make_unique()
    adata.obs_names_make_unique('-') 
    return adata

def preprocess_and_merge_adatas(adatas):

        # Step 2: Preprocess each AnnData object (example steps)
    for adata in adatas:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # Merge the AnnData objects into one
    adata_combined = adatas[0].concatenate(adatas[1:], batch_key='batch')
    return adata_combined




# Directly read the barcodes from the .tsv file and subset the AnnData object
def subset_anndata_from_barcodes_file(adata, barcodes_tsv_path):
    # Read barcodes into a pandas Series (squeeze=True converts the DataFrame to a Series)
    barcodes = pd.read_csv(barcodes_tsv_path, header=None, squeeze=True)

    # Subset the AnnData object
    return adata[adata.obs.index.isin(barcodes)].copy()




def subset_anndata_from_cluster_dictionary(adata, cluster_json_path):
    # Load the cluster dictionary from the JSON file
    with open(cluster_json_path) as json_file:
        cluster_dict = json.load(json_file)
    
    # Extract the barcodes (keys of the dictionary)
    barcodes = list(cluster_dict.keys())
    
    # Subset the AnnData object to only include the cells with the barcodes from the cluster dictionary
    adata_subset = adata[adata.obs.index.isin(barcodes)].copy()
    
    return adata_subset



def preprocess_adata(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    return adata



