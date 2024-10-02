using Serialization
using CSV
using DataFrames
using Glob
using Base.Threads

const OUTPUT_DIRECTORY = "/pool001/dswartz/viral_coev_mutation_migration_sweep_2"
const MIGRATION_RATES = exp10.(LinRange(-7, -0.5, 10))
const MUTATION_RATES = LinRange(0.001,0.02,10)

function calculate_survival_probability(migration_rate_idx, mutation_rate_idx)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$(migration_rate_idx)_mutation_rate_idx_$(mutation_rate_idx)")
    files = glob("*.jld2", output_subdirectory)
    num_replicates = length(files)
    survival_results = zeros(Bool, num_replicates)
    
    println((mutation_rate_idx, migration_rate_idx))
    flush(stdout)

    @threads for i in 1:num_replicates
        data = open(deserialize, files[i])
        total_infected_per_deme = data[1] # Assuming this is the correct structure
        survival_results[i] = sum(total_infected_per_deme[:, end]) > 0
    end

    return count(survival_results) / num_replicates
end

# DataFrame to store the results
df = DataFrame(MigrationRateIdx = Int[], MutationRateIdx = Int[], SurvivalProbability = Float64[])

for migration_rate_idx in 1:length(MIGRATION_RATES)
    for mutation_rate_idx in 1:length(MUTATION_RATES)
        prob = calculate_survival_probability(migration_rate_idx, mutation_rate_idx)
        push!(df, (migration_rate_idx, mutation_rate_idx, prob))
    end
end

# Save results to CSV file
csv_file = joinpath(OUTPUT_DIRECTORY, "survival_probabilities.csv")
CSV.write(csv_file, df)
