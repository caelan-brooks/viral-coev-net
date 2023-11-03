#!/bin/bash

#SBATCH -N 1                  # Number of nodes
#SBATCH -n 8                  # Number of cores
#SBATCH -p sched_mit_hill   # Partition name
#SBATCH --mem-per-cpu=3000    # Memory per CPU (in MB), adjust based on your script's requirements
#SBATCH -o analysis_output_%j.out # Output log (%j is replaced by job ID)
#SBATCH -e analysis_error_%j.err  # Error log

# Load Julia module and activate the project (adjust the paths as needed)
module load julia
export JULIA_PROJECT="/home/dswartz/viral-coev-net/julia_code/cluster_code/cluster_project"

# Execute your Julia script using the specific project and all available cores
julia --project=$JULIA_PROJECT -t auto preprocess_beta_sweep.jl
