using Serialization
using Glob
using Printf
using Statistics
using Base.Threads
include("/home/dswartz/viral-coev-net/julia_code/coevolution_network_base.jl")
using .CoevolutionNetworkBase

const OUTPUT_DIRECTORY = "/pool001/dswartz/migration_network_sweep_results"
const MIGRATION_RATES = exp10.(LinRange(-6, 0.5, 10))
const NETWORK_SIZES = 2:6
const LOCAL_RESULTS_DIRECTORY = "./migration_network_sweep_results"

println("Number of threads: ", nthreads())

function calculate_and_save_total_infected(network_size, migration_rate_idx)
    # Directory to save the total infected data
    save_directory = joinpath(LOCAL_RESULTS_DIRECTORY, "total_infected")
    if !isdir(save_directory)
        mkpath(save_directory)
    end
    
    # Define the directory and pattern for replicates
    dir = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(migration_rate_idx)")
    pattern = "replicate_*.jld2"

    # Check the number of replicates by counting the files
    files = glob(pattern, dir)
    num_replicates = length(files)

    # Initialize a collection to hold all the data
    all_replicates_data = [Dict() for _ in 1:num_replicates]

    @threads for replicate_idx in 1:num_replicates
        file_path = joinpath(dir, "replicate_$(replicate_idx).jld2")
        
        if isfile(file_path)
            try
                simulation = open(deserialize, file_path)
                total_infected_per_deme = calculate_total_infected_per_deme(simulation)
                times = simulation.duration_times  # Assuming time steps are in discrete units starting from 0
                
                # Insert the data for this replicate into the collection
                all_replicates_data[replicate_idx] = Dict("total_infected" => total_infected_per_deme, "times" => times)
                
            catch e
                println("Error processing file $(file_path): $e")
            end
        end
    end

    # Save all replicates data to a single file
    save_file_path = joinpath(save_directory, "network_size_$(network_size)_migration_rate_idx_$(migration_rate_idx).jld2")
    open(save_file_path, "w") do file
        serialize(file, all_replicates_data)
    end
end

function main(job_id_arg)
    job_id = parse(Int, job_id_arg)
    total_combinations = length(NETWORK_SIZES) * length(MIGRATION_RATES)
    
    # Check if job_id is within the valid range
    if job_id < 1 || job_id > total_combinations
        error("Invalid job ID: $job_id")
    end

    # Calculate indices for network size and migration rate
    network_size_idx = (job_id - 1) ÷ length(MIGRATION_RATES) + 1
    migration_rate_idx = (job_id - 1) % length(MIGRATION_RATES) + 1
    network_size = NETWORK_SIZES[network_size_idx]
    
    println("Processing for Network Size: $network_size, Migration Rate Index: $migration_rate_idx")
    calculate_and_save_total_infected(network_size, migration_rate_idx)
end

# Get job_id from command line arguments or environment variable
job_id_arg = get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
main(job_id_arg)
