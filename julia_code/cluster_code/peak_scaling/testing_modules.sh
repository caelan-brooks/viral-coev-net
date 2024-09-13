#!/bin/bash
#SBATCH -N 1                # Number of nodes
#SBATCH -n 1               # Number of cores
#SBATCH -p sched_mit_hill   # Partition name
#SBATCH --mem-per-cpu=10  # Memory per CPU (adjust this based on your requirements)
#SBATCH -t 0-00:10

# Path to the Julia project/environment
JULIA_PROJECT_PATH="/home/dswartz/viral-coev-net/julia_code/cluster_code/cluster_project"

source /etc/profile.d/modules.sh   # Adjust this if the file location is different on your system

# printenv

# echo $PATH

type module

# which module

module load julia

# Execute your Julia script using the specific project and all available cores
julia --project=$JULIA_PROJECT_PATH -e 'println("hello")'