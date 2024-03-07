using Serialization
using CSV
using DataFrames
using Glob
using Base.Threads
using Statistics

const OUTPUT_DIRECTORY = "/pool001/dswartz/real_network"
const CSV_OUTPUT_DIRECTORY = "/pool001/dswartz/real_network/csv_outputs"  # Directory for CSV outputs
const OUTBREAK_DEMES = collect(1:20)

println("Number of threads: ", nthreads())

# Create CSV output directory if it doesn't exist
isdir(CSV_OUTPUT_DIRECTORY) || mkdir(CSV_OUTPUT_DIRECTORY)

# Define read_data as a global function
function read_data(file)
    return open(deserialize, file)
end

# Initialize survival_probabilities as an array with NaN
survival_probabilities = fill(NaN, length(OUTBREAK_DEMES))

# Define the function to process each replicate
function process_replicate(file_path)
    try
        total_infected_per_deme, _ = read_data(file_path)
        total_infected = vec(sum(total_infected_per_deme, dims=1))

        survived = total_infected[end] > 0
        survived_flag = survived ? 1 : 0
        extinction_time = !survived ? findfirst(total_infected .== 0) : NaN
        
        return (survived_flag, extinction_time)
    catch e
        println("error reading file: error is $e")
        return (0, NaN) # Assuming survival flag as 0 and extinction time as NaN in case of error
    end
end

# Loop over each outbreak deme
for (idx, outbreak_deme) in enumerate(OUTBREAK_DEMES)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "outbreak_deme_idx_$idx")
    replicate_files = glob("*.jld2", output_subdirectory)
    num_replicates = length(replicate_files)

    println(num_replicates)
    flush(stdout)
    
    results = Array{Any}(undef, num_replicates)
    
    # Parallel processing of each replicate
    @threads for file_idx in 1:num_replicates
        results[file_idx] = process_replicate(replicate_files[file_idx])
    end

    survived_results = [result[1] for result in results]
    extinction_times = [result[2] for result in results]
    
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




