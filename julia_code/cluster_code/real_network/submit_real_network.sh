#!/bin/bash
#SBATCH -N 1                # Number of nodes
#SBATCH -n 6               # Number of cores
#SBATCH -p sched_mit_hill   # Partition name
#SBATCH -t 0-11:58          # time taken is 4 hrs
#SBATCH --mem-per-cpu=2000  # Memory per CPU (adjust this based on your requirements)
#SBATCH --array=1-20        # Create a job array for 10 jobs

# Path to the Julia project/environment
JULIA_PROJECT_PATH="/home/dswartz/viral-coev-net/julia_code/cluster_code/cluster_project"

module add /home/software/modulefiles/julia/

# Execute your Julia script using the specific project and all available cores
julia --project=$JULIA_PROJECT_PATH -t 12 real_network_sweep.jl