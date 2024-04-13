#!/bin/bash

TEMPLATE_FILE="template_structure.yaml"  # The path must be correct


RUN_ID="$1"
SAMPLE_ID="$2" #I believe this takes param

RUNS_DIR_ID="./runs/${RUN_ID}"

# Function to create directories based on the template
create_directories() {
    local run_id=$1
    local dirs

    # Extract the keys as directory names using yq from the placeholder
    dirs=$(yq e '.run_id_placeholder | keys' "$TEMPLATE_FILE" | awk '{print $2}')  # Extract keys without dashes

    # Create the run directory
    mkdir -p "${RUNS_DIR_ID}/${SAMPLE_ID}"

    # Create subdirectories
    for dir in $dirs; do
        # Remove any potential surrounding quotes
        dir=$(echo $dir | tr -d '"')
        # Check if dir is not empty and not a comment
        if [[ -n "$dir" && "$dir" != \#* ]]; then
            mkdir -p "${RUNS_DIR_ID}/${SAMPLE_ID}/${dir}"
        fi
    done

    echo "Directory structure initialized for run: ${run_id}"
}

# Ensure the TEMPLATE_FILE exists
if [ ! -f "$TEMPLATE_FILE" ]; then
    echo "The template structure file does not exist: $TEMPLATE_FILE"
    exit 1
fi

# Ensure a run ID is provided
if [ -z "$SAMPLE_ID" ]; then
    echo "Usage: $0 SAMPLE_ID"
    exit 1
fi

create_directories "${RUNS_DIR_ID}${SAMPLE_ID}"
