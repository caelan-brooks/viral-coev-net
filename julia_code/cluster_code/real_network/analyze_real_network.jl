using Serialization
using CSV
using DataFrames
using Glob
using Base.Threads
using Statistics

const OUTPUT_DIRECTORY = "/pool001/dswartz/real_network"
const CSV_OUTPUT_DIRECTORY = "/pool001/dswartz/real_network/csv_outputs"  # Directory for CSV outputs
const OUTBREAK_DEMES = collect(1:20)

# Create CSV output directory if it doesn't exist
isdir(CSV_OUTPUT_DIRECTORY) || mkdir(CSV_OUTPUT_DIRECTORY)

# Define read_data as a global function
function read_data(file)
    return open(deserialize, file)
end

# Initialize survival_probabilities as an array with NaN
survival_probabilities = fill(NaN, length(OUTBREAK_DEMES))

# Loop over each migration rate
for (idx, outbreak_deme) in enumerate(OUTBREAK_DEMES)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "outbreak_deme_idx_$idx")
    replicate_files = glob("*.jld2", output_subdirectory)
    num_replicates = length(replicate_files)

    println(num_replicates)
    flush(stdout)
    
    survived_results = zeros(Int, num_replicates)
    extinction_times = fill(NaN, num_replicates)
    
    # Process each replicate
    @threads for file_idx = 1:num_replicates
        try
            total_infected_per_deme, _ = read_data(replicate_files[file_idx])
            
            total_infected = vec(sum(total_infected_per_deme, dims=1))
            
            # Check if the pathogen survived
            survived = total_infected[end] > 0
            survived_results[file_idx] = survived ? 1 : 0
            if !survived
                extinction_times[file_idx] = findfirst(total_infected .== 0)
            end
        catch e
            println("error reading file: error is $e")
        end
    end
    
    valid_extinction_times = extinction_times[.!isnan.(extinction_times)]
    
    println("Maximum Extinction Time: ", maximum(valid_extinction_times))
    println("Average Extinction Time: ", mean(valid_extinction_times))
    println("number of times greater than 60: ", count(valid_extinction_times .> 60))
    println("number times greater than 70: ", count(valid_extinction_times .> 70))
    println("number times greater than 80: ", count(valid_extinction_times .> 80))
    
    # Update survival probability for this migration rate in the array
    survival_probabilities[idx] = sum(survived_results) / num_replicates
end


# Convert to DataFrame for CSV
df = DataFrame(
    OutbreakDeme = OUTBREAK_DEMES,
    SurvivalProbability = survival_probabilities
)

# Save to CSV
CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "analysis_results.csv"), df)




