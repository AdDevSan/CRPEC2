import os
import sys
import yaml

def create_directories(run_id, sample_id, template_file):
    # Load the YAML file
    try:
        with open(template_file, 'r') as file:
            data = yaml.safe_load(file)
    except Exception as e:
        print(f"Error reading the YAML file: {e}")
        sys.exit(1)

    # Extract keys from the 'run_id_placeholder' section
    directories = data['run_id_placeholder'].keys()

    # Base path for directories
    base_path = f"./runs/{run_id}/{sample_id}"

    # Create the run and sample directories
    os.makedirs(base_path, exist_ok=True)

    # Create subdirectories
    for directory in directories:
        os.makedirs(os.path.join(base_path, directory), exist_ok=True)

    print(f"Directory structure initialized for run: {run_id}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python create_directories.py <run_id> <sample_id> <template_file>")
        sys.exit(1)

    run_id = sys.argv[1]
    sample_id = sys.argv[2]
    template_file = sys.argv[3]

    create_directories(run_id, sample_id, template_file)

