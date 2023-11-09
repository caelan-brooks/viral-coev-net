#!/bin/bash
#SBATCH -N 1                # Number of nodes
#SBATCH -n 1                # Number of cores
#SBATCH -p sched_mit_hill   # Partition name
#SBATCH --mem=1000          # Adjust memory based on your requirements
#SBATCH --array=1-50        # 5 network sizes times 10 migration rates

# Echo the current working directory
echo "Current working directory: $(pwd)"

# Base directory
BASE_DIR="/pool001/dswartz/migration_network_sweep_results/"

# Calculate network size and migration rate index based on SLURM_ARRAY_TASK_ID
NETWORK_SIZE=$((2 + (SLURM_ARRAY_TASK_ID - 1) / 10))
MIGRATION_IDX=$(((SLURM_ARRAY_TASK_ID - 1) % 10 + 1))

# Define the directory to be emptied
DIR_TO_EMPTY="$BASE_DIR/network_size_${NETWORK_SIZE}/migration_rate_idx_${MIGRATION_IDX}"

# Path to the empty directory
EMPTY_DIR="$(pwd)/empty_directory"

# Use rsync to delete all files in the target directory
rsync --recursive --delete --verbose "$EMPTY_DIR/" "$DIR_TO_EMPTY"

echo "Contents of $DIR_TO_EMPTY have been deleted."
