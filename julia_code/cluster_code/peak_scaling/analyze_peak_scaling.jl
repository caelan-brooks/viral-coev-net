using Serialization
using CSV
using DataFrames
using Glob
using Base.Threads
using Statistics

const OUTPUT_DIRECTORY = "/pool001/dswartz/peak_scaling"

df = CSV.read(joinpath(OUTPUT_DIRECTORY,"migration_rates.csv"), DataFrame)
const MIGRATION_RATES = df.MIGRATION_RATES

df = CSV.read(joinpath(OUTPUT_DIRECTORY, "host_population_sizes.csv"), DataFrame)
const HOST_POPULATION_SIZES = df.HOST_POPULATION_SIZES

println("Number of threads: ", nthreads())

function readata(file)
    return open(deserialize, file)
end


function calculate_survival_probability(migration_rate_idx, host_per_deme_idx)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "host_per_deme_idx_$(host_per_deme_idx)", "migration_rate_idx_$(migration_rate_idx)")
    files = glob("*.jld2", output_subdirectory)
    num_replicates = length(files)
    survival_results = zeros(Bool, num_replicates)
    
    println((host_per_deme_idx, migration_rate_idx, num_replicates))
    flush(stdout)

    @threads for i in 1:num_replicates
        total_infected_per_deme, _ = open(deserialize, files[i])
        # total_infected_per_deme = data[1] # Assuming this is the correct structure
        survival_results[i] = sum(total_infected_per_deme[:, end]) > 0
    end

    return count(survival_results) / num_replicates
end

df = DataFrame(MigrationRateIdx = Int[], HostPerDemeIdx = Int[], SurvivalProbability = Float64[])

for host_idx = 1:length(HOST_POPULATION_SIZES)
    for migration_idx = 1:length(MIGRATION_RATES)
        prob = calculate_survival_probability(migration_idx, host_idx)
        push!(df, (migration_idx, host_idx, prob))
    end
end

filepath = joinpath(OUTPUT_DIRECTORY, "survival_probabilities.csv")
CSV.write(filepath, df)