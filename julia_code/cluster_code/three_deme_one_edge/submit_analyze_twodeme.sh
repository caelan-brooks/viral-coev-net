#!/bin/bash
#SBATCH -N 1                # Number of nodes
#SBATCH -n 1               # Number of cores
#SBATCH -p sched_mit_hill   # Partition name
#SBATCH --mem-per-cpu=1000  # Memory per CPU (adjust this based on your requirements)

# Path to the Julia project/environment
JULIA_PROJECT_PATH="/home/dswartz/viral-coev-net/julia_code/cluster_code/cluster_project"

module add /home/software/modulefiles/julia/

# Execute your Julia script using the specific project and all available cores
julia --project=$JULIA_PROJECT_PATH -t auto analyze_twodeme.jl
