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

function analyze_data(network_size, migration_rate_idx, num_replicates)
    total_survivals = 0
    
    for replicate_idx in 1:num_replicates
        file_path = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(migration_rate_idx)", "replicate_$(replicate_idx).jld2")
        
        if isfile(file_path)
            try
                simulation = open(deserialize, file_path)
                total_infected = calculate_total_infected(simulation)
                if total_infected[end] > 0
                    total_survivals += 1
                end
            catch e
                println("Error processing file $(file_path): $e")
            end
        end
    end

    survival_probability = total_survivals / num_replicates
    return survival_probability, total_survivals
end

# Preallocate arrays for results
survival_probabilities = zeros(length(NETWORK_SIZES), length(MIGRATION_RATES))
num_replicates = zeros(Int, length(NETWORK_SIZES), length(MIGRATION_RATES))

# Loop over network sizes and migration rates
for (size_idx, network_size) in enumerate(NETWORK_SIZES)
    for (rate_idx, migration_rate) in enumerate(MIGRATION_RATES)
        # Define the directory and pattern
        dir = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(rate_idx)")
        pattern = "replicate_*.jld2"

        # Check if the directory exists
        if isdir(dir)
            # Detect the number of replicates by counting the files
            files = glob(pattern, dir)
            num_replicates[size_idx, rate_idx] = length(files)
            println("Network Size: $network_size, Migration Rate Index: $rate_idx, Replicates: $(num_replicates[size_idx, rate_idx])")
            flush(stdout)
            
            # Analyze data
            survival_probabilities[size_idx, rate_idx], total_survivals = analyze_data(network_size, rate_idx, num_replicates[size_idx, rate_idx])
            
            println("Total Survivals: $total_survivals, Survival Probability: $(survival_probabilities[size_idx, rate_idx])")
        else
            println("Directory does not exist: $dir")
        end
    end
end

# Save the results
if !isdir(LOCAL_RESULTS_DIRECTORY)
    mkdir(LOCAL_RESULTS_DIRECTORY)
end
save_path = joinpath(LOCAL_RESULTS_DIRECTORY, "migration_network_sweep_results.jld2")
open(save_path, "w") do file
    serialize(file, Dict("survival_probabilities" => survival_probabilities, "num_replicates" => num_replicates))
end
