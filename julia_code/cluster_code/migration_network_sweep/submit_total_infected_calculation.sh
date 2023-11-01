#!/bin/bash
#SBATCH -N 1                # Number of nodes
#SBATCH -n 8               # Number of cores
#SBATCH -p sched_mit_hill   # Partition name
#SBATCH --mem-per-cpu=4000  # Memory per CPU (adjust this based on your requirements)
#SBATCH --array=1-50        # 5 network sizes times 10 migration rates

# Path to the Julia project/environment
JULIA_PROJECT_PATH="/home/dswartz/viral-coev-net/julia_code/cluster_code/cluster_project"

module load julia

# Execute your Julia script using the specific project and all available cores
julia --project=$JULIA_PROJECT_PATH -t auto analyze_infected_trajectories.jl
