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
    # Preallocate an array to store the survival status for each replicate
    survival_statuses = zeros(Bool, num_replicates)

    @threads for replicate_idx in 1:num_replicates
        file_path = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(migration_rate_idx)", "replicate_$(replicate_idx).jld2")

        if isfile(file_path)
            try
                simulation = open(deserialize, file_path)
                total_infected = calculate_total_infected(simulation)
                survival_statuses[replicate_idx] = total_infected[end] > 0
            catch e
                println("Error processing file $(file_path): $e")
            end
        end
    end

    total_survivals = sum(survival_statuses)
    survival_probability = total_survivals / num_replicates
    return survival_probability, total_survivals
end

# # Preallocate arrays for results
# survival_probabilities = zeros(length(NETWORK_SIZES), length(MIGRATION_RATES))
# num_replicates = zeros(Int, length(NETWORK_SIZES), length(MIGRATION_RATES))

# # Loop over network sizes and migration rates
# for (size_idx, network_size) in enumerate(NETWORK_SIZES)
#     for (rate_idx, migration_rate) in enumerate(MIGRATION_RATES)
#         # Define the directory and pattern
#         dir = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(rate_idx)")
#         pattern = "replicate_*.jld2"

#         # Check if the directory exists
#         if isdir(dir)
#             # Detect the number of replicates by counting the files
#             files = glob(pattern, dir)
#             num_replicates[size_idx, rate_idx] = length(files)
#             println("Network Size: $network_size, Migration Rate Index: $rate_idx, Replicates: $(num_replicates[size_idx, rate_idx])")
#             flush(stdout)
            
#             # Analyze data
#             survival_probabilities[size_idx, rate_idx], total_survivals = analyze_data(network_size, rate_idx, num_replicates[size_idx, rate_idx])
            
#             println("Total Survivals: $total_survivals, Survival Probability: $(survival_probabilities[size_idx, rate_idx])")
#         else
#             println("Directory does not exist: $dir")
#         end
#     end
# end

# # Save the results
# if !isdir(LOCAL_RESULTS_DIRECTORY)
#     mkdir(LOCAL_RESULTS_DIRECTORY)
# end
# save_path = joinpath(LOCAL_RESULTS_DIRECTORY, "migration_network_sweep_results.jld2")
# open(save_path, "w") do file
#     serialize(file, Dict("survival_probabilities" => survival_probabilities, "num_replicates" => num_replicates))
# end

# Function to calculate and save total infected individuals
function calculate_and_save_total_infected(network_size, migration_rate_idx, num_replicates)
    # Directory to save the total infected data
    save_directory = joinpath(LOCAL_RESULTS_DIRECTORY, "total_infected", "network_size_$(network_size)", "migration_rate_idx_$(migration_rate_idx)")
    if !isdir(save_directory)
        mkpath(save_directory)
    end
    
    @threads for replicate_idx in 1:num_replicates
        file_path = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(migration_rate_idx)", "replicate_$(replicate_idx).jld2")
        
        if isfile(file_path)
            try
                simulation = open(deserialize, file_path)
                total_infected = calculate_total_infected(simulation)
                times = simulation.duration_times  # Assuming time steps are in discrete units starting from 0
                
                # Save the total infected data
                save_file_path = joinpath(save_directory, "replicate_$(replicate_idx).jld2")
                open(save_file_path, "w") do file
                    serialize(file, Dict("total_infected" => total_infected, "times" => times))
                end
                
                
            catch e
                println("Error processing file $(file_path): $e")
            end
        end
    end
end

# Calculate and save total infected individuals for each network size and migration rate
for (size_idx, network_size) in enumerate(NETWORK_SIZES)
    for (rate_idx, migration_rate) in enumerate(MIGRATION_RATES)
        dir = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(rate_idx)")
        
        if isdir(dir)
            files = glob("replicate_*.jld2", dir)
            num_replicates = length(files)
            println("Network Size: $network_size, Migration Rate Index: $rate_idx, Replicates: $num_replicates")
            flush(stdout)
            
            calculate_and_save_total_infected(network_size, rate_idx, num_replicates)
        else
            println("Directory does not exist: $dir")
        end
    end
end
