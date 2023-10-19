using Serialization
using Glob
using Printf
using Statistics

include("../coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DATA_DIRECTORY = "/pool001/dswartz/viral_coev_output/betasweep_results"
const PROCESSED_DATA_DIRECTORY = "./processed_data_betasweep"
const BETA_VALUES = LinRange(1.05, 6, 10)  # Example beta values to sweep over

# Create the processed data directory if it doesn't exist
if !isdir(PROCESSED_DATA_DIRECTORY)
    mkdir(PROCESSED_DATA_DIRECTORY)
end

"""
    classify_extinction_type(simulation)

Determine the type of extinction in the simulation.
Returns:
- "demographic" for extinction due to demographic noise
- "immune" for extinction due to immune response
- "survival" if no extinction
- "unknown" in case of errors
"""
function classify_extinction_type(simulation)
    try
        total_infected_at_t4 = calculate_total_infected(simulation)[Int(round(4 / simulation.dt))]
        total_infected_at_end = calculate_total_infected(simulation)[end]

        if total_infected_at_t4 == 0
            return "demographic"
        elseif total_infected_at_end == 0
            return "immune"
        else
            return "survival"
        end
    catch e
        println("Error processing data: $e")
        return "unknown"
    end
end

"""
    analyze_data_for_beta(beta_value)

Analyze the data for a given beta value and return counts of each extinction type.
"""
function analyze_data_for_beta(beta_value)
    files = glob("simulation_results_beta_$(beta_value)_replicate_*.jld2", DATA_DIRECTORY)
    outcomes = Dict("demographic" => 0, "immune" => 0, "survival" => 0, "unknown" => 0)

    Threads.@threads for idx in 1:length(files)
        file = files[idx]
        try
            simulation = open(deserialize, file)
            extinction_type = classify_extinction_type(simulation)
            outcomes[extinction_type] += 1
        catch e
            println("Error processing file $file: $e")
            outcomes["unknown"] += 1
        end
    end

    return outcomes
end

# Analyzing the data and saving the results
results = Dict()
for (idx, beta) in enumerate(BETA_VALUES)
    println(beta)
    results[idx] = analyze_data_for_beta(beta)
end


# Saving the results
save_path = joinpath(PROCESSED_DATA_DIRECTORY, "extinction_results.jld2")
open(save_path, "w") do file
    serialize(file, results)
end
