using Serialization
using CSV
using DataFrames
using Glob

const OUTPUT_DIRECTORY = "/pool001/dswartz/real_network"
const CSV_OUTPUT_DIRECTORY = "/pool001/dswartz/real_network/csv_outputs"  # Directory for CSV outputs
const OUTBREAK_DEMES = collect(1:20)

# Create CSV output directory if it doesn't exist
isdir(CSV_OUTPUT_DIRECTORY) || mkdir(CSV_OUTPUT_DIRECTORY)

# Define read_data as a global function
function read_data(file)
    open(deserialize, file)
end

# Initialize survival_probabilities as an array with NaN
survival_probabilities = fill(NaN, length(OUTBREAK_DEMES))

# Loop over each migration rate
for (idx, outbreak_deme) in enumerate(OUTBREAK_DEMES)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "outbreak_deme_idx_$idx")
    replicate_files = glob("*.jld2", output_subdirectory)
    num_replicates = length(replicate_files)
    num_survived = 0
    println(num_replicates)
    flush(stdout)

    # Process each replicate
    for file in replicate_files
        total_infected_per_deme, _ = read_data(file)
        
        # Check if the pathogen survived
        survived = sum(total_infected_per_deme[:, end]) > 0
        num_survived += survived ? 1 : 0
    end

    # Update survival probability for this migration rate in the array
    survival_probabilities[idx] = num_survived / num_replicates
end


# Convert to DataFrame for CSV
df = DataFrame(
    OutbreakDeme = OUTBREAK_DEMES,
    SurvivalProbability = survival_probabilities
)

# Save to CSV
CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "analysis_results.csv"), df)




