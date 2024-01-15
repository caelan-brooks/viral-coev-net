using Serialization
using CSV
using DataFrames
using Glob
using Base.Threads

const OUTPUT_DIRECTORY = "/pool001/dswartz/twodeme_final"
const MIGRATION_RATES = [0; exp10.(LinRange(-7, -0.5, 10)); 0]

# Initialize survival_probabilities as an array with NaN
survival_probabilities = fill(NaN, length(MIGRATION_RATES))

# Loop over each migration rate
for (idx, migration_rate) in enumerate(MIGRATION_RATES)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$idx")
    replicate_files = glob("*.jld2", output_subdirectory)  # Assuming files are saved with .jld2 extension
    num_replicates = length(replicate_files)
    num_survived = 0

    # Process each replicate
    @threads for file in replicate_files
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
    MigrationRate = MIGRATION_RATES,
    SurvivalProbability = survival_probabilities
)

# Save to CSV
CSV.write("analysis_results.csv", df)
