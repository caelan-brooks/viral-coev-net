using Serialization
using CSV
using DataFrames
using Glob
using Base.Threads
using Statistics

const OUTPUT_DIRECTORY = "/pool001/dswartz/pop_size_one_deme"
const CSV_OUTPUT_DIRECTORY = joinpath(OUTPUT_DIRECTORY, "csv_outputs")  # Directory for CSV outputs
const HOST_POPULATION_SIZES = exp10.(LinRange(5.0,8.0,10))

# Create CSV output directory if it doesn't exist
isdir(CSV_OUTPUT_DIRECTORY) || mkdir(CSV_OUTPUT_DIRECTORY)

println("Number of threads: ", nthreads())

# Define read_data as a global function
function read_data(file)
    open(deserialize, file)
end

# Function to process each replicate
function process_replicate(file_path::String)
    try
        total_infected_per_deme, duration_times = read_data(file_path)
        total_infected = vec(total_infected_per_deme) # only one deme

        maximum_infected_deme_1 = maximum(total_infected_per_deme[1,:])

        survived = total_infected[end] > 0
        survived_flag = survived ? 1 : 0
        extinction_time = !survived ? duration_times[findfirst(total_infected .== 0)] : NaN
        
        return (survived_flag, extinction_time, maximum_infected_deme_1)
    catch e
        println("error reading file: error is $e")
        return (0, NaN, NaN, NaN) # Assuming default values in case of error
    end
end

function analyze_all_results(HOST_POPULATION_SIZES, OUTPUT_DIRECTORY::String)
    survival_probabilities = fill(NaN, length(HOST_POPULATION_SIZES))

    # Your existing process_replicate function remains unchanged

    for idx in eachindex(HOST_POPULATION_SIZES)
        output_subdirectory = joinpath(OUTPUT_DIRECTORY, "host_per_deme_idx_$idx")
        replicate_files = glob("*.jld2", output_subdirectory)
        num_replicates = length(replicate_files)
        println(num_replicates)
        flush(stdout)

        results = Array{Any}(undef, num_replicates)

        @time begin
            @threads for file_idx in 1:num_replicates
                results[file_idx] = process_replicate(replicate_files[file_idx])
            end
        end

        survived_results = [result[1] for result in results]
        extinction_times = [result[2] for result in results]
        maximum_infected_deme_1 = [result[3] for result in results]

        valid_extinction_times = extinction_times[.!isnan.(extinction_times)]

        println("Maximum Extinction Time: ", maximum(valid_extinction_times))
        println("Average Extinction Time: ", mean(valid_extinction_times))
        println("number of times greater than 60: ", count(valid_extinction_times .> 60))
        println("number times greater than 70: ", count(valid_extinction_times .> 70))
        println("number times greater than 80: ", count(valid_extinction_times .> 80))
        println("number times greater than 110: ", count(valid_extinction_times .> 110))
        println("Peak infected in deme 1: ", mean(maximum_infected_deme_1), " +- ", std(maximum_infected_deme_1))

        survival_probabilities[idx] = sum(survived_results) / num_replicates
    end

    return survival_probabilities
end

survival_probabilities = analyze_all_results(HOST_POPULATION_SIZES, OUTPUT_DIRECTORY)

# Convert to DataFrame for CSV
df = DataFrame(
    HostPopulationSize = HOST_POPULATION_SIZES,
    SurvivalProbability = survival_probabilities
)

# Save to CSV
CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "analysis_results.csv"), df)

