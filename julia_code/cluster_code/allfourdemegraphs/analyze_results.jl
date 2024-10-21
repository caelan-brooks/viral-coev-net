using Serialization
using Glob
using Printf
using Statistics
using Base.Threads
include("/home/dswartz/viral-coev-net/julia_code/coevolution_network_base.jl")
using .CoevolutionNetworkBase

include("all_adjacency_matrices.jl")

const DATA_DIRECTORY = "/pool001/dswartz/all_four_deme_unnorm_PL_with_dx"
const LOCAL_RESULTS_DIRECTORY = "./plotting_data"

const MIGRATION_RATES = exp10.(LinRange(-6, -0.5, 10))

println("Number of threads: ", nthreads())

function analyze_data(adjacency_matrix_idx, migration_rate_idx, num_replicates; cutoff=0)
    # Preallocate an array to store the survival status for each replicate
    survival_statuses = zeros(Bool, num_replicates)

    @threads for replicate_idx in 1:num_replicates
        file_path = joinpath(DATA_DIRECTORY, "adjacency_matrix_idx_$(adjacency_matrix_idx)", "migration_rate_idx_$(migration_rate_idx)", "replicate_$(replicate_idx).jld2")

        if isfile(file_path)
            try
                total_infected_per_deme = open(deserialize,file_path)
                total_infected_end = sum(total_infected_per_deme[:, end]) 
                survival_statuses[replicate_idx] = total_infected_end > cutoff
            catch e
                println("Error processing file $(file_path): $e")
            end
        end
    end

    total_survivals = sum(survival_statuses)
    survival_probability = total_survivals / num_replicates
    return survival_probability, total_survivals
end

# Initialize a dictionary to store the results
results = Dict()

# Iterate over each adjacency matrix and migration rate
for adjacency_matrix_idx in 1:length(all_adjacency_matrices)
    for (rate_idx, migration_rate) in enumerate(MIGRATION_RATES)
        dir = joinpath(DATA_DIRECTORY, "adjacency_matrix_idx_$(adjacency_matrix_idx)", "migration_rate_idx_$(rate_idx)")
        
        if isdir(dir)
            files = glob("replicate_*.jld2", dir)
            num_replicates = length(files)
            println("Adjacency Matrix: $adjacency_matrix_idx, Migration Rate Index: $rate_idx, Replicates: $num_replicates")
            flush(stdout)
            
            # Call the analysis function and store the results
            survival_probability, total_survivals = analyze_data(adjacency_matrix_idx, rate_idx, num_replicates)
            println(survival_probability)
            results[(adjacency_matrix_idx, rate_idx)] = (survival_probability, total_survivals, num_replicates)
        else
            println("Directory does not exist: $dir")
        end
    end
end

# Serialize the results dictionary to a file
output_file = joinpath(LOCAL_RESULTS_DIRECTORY, "analysis_results.jld2")
open(output_file, "w") do file
    serialize(file, results)
end

println("Analysis results saved to $output_file")
