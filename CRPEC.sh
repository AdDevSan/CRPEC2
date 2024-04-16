#!/bin/bash
#SBATCH --job-name=CRPEC
#SBATCH --partition=gpu3090
#SBATCH --nodes=1
#SBATCH --qos=gpu3090
#SBATCH --cpus-per-task=8
#SBATCH --output=CRPEC%j.out
#SBATCH --error=CRPEC%j.err



# Generate a run ID based on the current date and time


# Directory where trident directories are stored
INPUT_DIR="input"
RUNS_DIR="./runs"
# Generate a unique run ID
RUN_ID="CRPEC_run_$(date +%Y%m%d_%H%M%S)"

RUNS_DIR_ID="${RUNS_DIR}/${RUN_ID}"
mkdir "${RUNS_DIR_ID}"  # runs/CRPEC_run_X

# Loop over each trident directory in the input directory
for TRIDENT_DIRECTORY in ${INPUT_DIR}/*; do
    # Extract the base name for the trident directory to use as SAMPLE_ID
    SAMPLE_ID=$(basename "${TRIDENT_DIRECTORY}")
    

    

    # Initialize directories for the new run
    bash initialize_directories.sh "${RUN_ID}" "${SAMPLE_ID}"

    # Define directory paths
    SAMPLE_DIR="${RUNS_DIR_ID}/${SAMPLE_ID}"
    INITIAL_CLUSTERS_DIR="${RUNS_DIR_ID}/${SAMPLE_ID}/initial_clusters"
    REFINED_CLUSTERS_DIR="${RUNS_DIR_ID}/${SAMPLE_ID}/refined_clusters"
    PREDICTED_CLUSTERS_DIR="${RUNS_DIR_ID}/${SAMPLE_ID}/predicted_clusters"
    SAMPLE_200_DIR="${RUNS_DIR_ID}/${SAMPLE_ID}/sample_200"





    LOOP_RUNS=100


    #python trident_preprocess_to_h5ad_barcodes.py -t <trident_directory> -hd <h5ad_directory> -bd <filtered_barcodes_directory>
    python3 trident_preprocess_to_h5ad_barcodes.py -t "${TRIDENT_DIRECTORY}" -hd "${SAMPLE_DIR}" -bd "${SAMPLE_DIR}"

    #OUTPUTS FROM TRIDENT PREPROCESS:
    H5AD_FILE_PATH_ORIGINAL="${RUNS_DIR_ID}/${SAMPLE_ID}/full_dataset_processed.h5ad"
    H5AD_FILE_PATH="${RUNS_DIR_ID}/${SAMPLE_ID}/full_dataset_1000var_pca_processed.h5ad"
    FILTERED_BARCODES="${SAMPLE_DIR}/filtered_barcodes.tsv"

    #'-b''-n''-s''-r'
    #args.barcodes_tsv, args.n_samples, args.output_dir_samples, args.output_dir_remainders
    python3 sample_200.py -b "${FILTERED_BARCODES}" -n ${LOOP_RUNS} -s "${SAMPLE_DIR}/sample_200" -r "${SAMPLE_DIR}/sample_remainder"



    Rscript R_scripts/h5ad_to_h5seurat.R --input "${H5AD_FILE_PATH_ORIGINAL}"
    #OUTPUT H5 SEURAT
    H5SEURAT_FILE_PATH="${RUNS_DIR_ID}/${SAMPLE_ID}/full_dataset_processed.h5seurat"


    Rscript preprocess_sc3.R --input.h5seurat.file.path "${H5SEURAT_FILE_PATH}" --input.barcodes.directory.path "${SAMPLE_200_DIR}" --output.directory "${INITIAL_CLUSTERS_DIR}"





    # REFINE CLUSTERS refine_clusters.sh
    bash refine_clusters.sh "${H5AD_FILE_PATH}" "${INITIAL_CLUSTERS_DIR}" "${REFINED_CLUSTERS_DIR}"

    # PREDICT CLUSTERS 
    bash predict_clusters.sh "${H5AD_FILE_PATH}" "${REFINED_CLUSTERS_DIR}" "${PREDICTED_CLUSTERS_DIR}"

    # GET ENSEMBLE CLUSTERS
    python3 get_consensus_final_cluster_dict.py -d "${SAMPLE_DIR}"

    ENSEMBLE_CLUSTER_PATH="${SAMPLE_DIR}/ensemble_cluster_full.json"

    python3 get_ensemble_adatas.py \
    --sample_name "${SAMPLE_ID}" \
    --ensemble_cluster_path "${ENSEMBLE_CLUSTER_PATH}" \
    --processed_adata_path "${H5AD_FILE_PATH_ORIGINAL}" \
    --output_directory_path "${SAMPLE_DIR}/ensemble_adatas"

    #python process_adata.py 'runs/CRPEC_run_20240414_031735' 'gene_profiles_mean.csv'
    

    python get_gene_profiles_dtf.py --run_id "${RUNS_DIR_ID}" --output_file "${RUNS_DIR_ID}/gene_profiles_mean.csv"
done
