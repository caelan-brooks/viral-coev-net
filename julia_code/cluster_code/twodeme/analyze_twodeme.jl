using Serialization
using CSV
using DataFrames
using Glob
using Base.Threads
using Statistics

const OUTPUT_DIRECTORY = "/pool001/dswartz/twodeme_convergent_small_dx"
const CSV_OUTPUT_DIRECTORY = "/pool001/dswartz/twodeme_convergent_small_dx/csv_outputs"  # Directory for CSV outputs
const MIGRATION_RATES = [0; exp10.(LinRange(-7, -0.5, 10)); 0]

# Create CSV output directory if it doesn't exist
isdir(CSV_OUTPUT_DIRECTORY) || mkdir(CSV_OUTPUT_DIRECTORY)

println("Number of threads: ", nthreads())

# Define read_data as a global function
function read_data(file)
    open(deserialize, file)
end

function process_and_save_histograms(migration_rate_idx::Int64; cutoff = 100)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$migration_rate_idx")
    replicate_files = glob("*.jld2", output_subdirectory)

    antigenic_variance_deme1 = []
    antigenic_variance_deme2 = []
    peak_time_difference = []
    variance_difference = []

    number_late_trajectories = 0

    for file in replicate_files
        total_infected_per_deme, antigenic_variance_per_deme = read_data(file)
        
        # Find indices of maximal infection for each deme because dt = 1/20 and thin by = 20 index corrresponds to "real" time!
        max_infected_idx_deme1 = argmax(total_infected_per_deme[1, :])
        min_infected_idx_deme1 = argmin(total_infected_per_deme[1, 6:end])
        
        max_infected_idx_deme2 = argmax(total_infected_per_deme[2, :])
        
        
        
        # Check if maximal infection is greater than cutoff before recording the antigenic variance
        if total_infected_per_deme[1, max_infected_idx_deme1] > cutoff
            push!(antigenic_variance_deme1, antigenic_variance_per_deme[1, max_infected_idx_deme1])
        end
        
        # if total_infected_per_deme[2, max_infected_idx_deme2] > cutoff
        #     push!(antigenic_variance_deme2, antigenic_variance_per_deme[2, max_infected_idx_deme2])
        #     push!(peak_time_difference, max_infected_idx_deme2 - max_infected_idx_deme1)
        #     push!(variance_difference, antigenic_variance_per_deme[2, max_infected_idx_deme2] - antigenic_variance_per_deme[1, max_infected_idx_deme1])
        # end
        
        first_cross_idx_deme2 = findfirst(total_infected_per_deme[2, :] .>= cutoff)
        
        if !isnothing(first_cross_idx_deme2)
            if first_cross_idx_deme2 > min_infected_idx_deme1
                number_late_trajectories += 1
            end
     
            if total_infected_per_deme[2, max_infected_idx_deme2] > cutoff && first_cross_idx_deme2 <= min_infected_idx_deme1
                push!(antigenic_variance_deme2, antigenic_variance_per_deme[2, max_infected_idx_deme2])
                push!(peak_time_difference, max_infected_idx_deme2 - max_infected_idx_deme1)
                push!(variance_difference, antigenic_variance_per_deme[2, max_infected_idx_deme2] - antigenic_variance_per_deme[1, max_infected_idx_deme1])
            end
        end
    end
    
    println("Number of late trajectories: ", number_late_trajectories)

    # Save the antigenic variance data
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "antigenic_variance_deme1_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(AntigenicVariance=antigenic_variance_deme1))
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "antigenic_variance_deme2_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(AntigenicVariance=antigenic_variance_deme2))
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "peak_time_difference_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(PeakTimeDifference=peak_time_difference))
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "variance_difference_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(VarianceDifference=variance_difference))
end

function process_trajectories(migration_rate_idx::Int64, num_replicates::Int64)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$migration_rate_idx")
    replicate_files = glob("*.jld2", output_subdirectory)
    
    # Initialize a matrix to store all trajectories
    trajectories_matrix = []

    for i in 1:min(num_replicates, length(replicate_files))
        file = replicate_files[i]
        total_infected_per_deme, _ = read_data(file)
        total_infected = sum(total_infected_per_deme, dims=1)  # Sum over the first index
        push!(trajectories_matrix, total_infected[:]')
    end

    # Convert the array of trajectories into a DataFrame
    traj_df = DataFrame(trajectories_matrix)

    # Save to CSV
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "trajectories_migration_rate_idx_$migration_rate_idx.csv"), traj_df)
end

# Function to process each replicate
function process_replicate(file_path::String)
    try
        total_infected_per_deme, _ = read_data(file_path)
        total_infected = vec(sum(total_infected_per_deme, dims=1))

        maximum_infected_deme_1 = maximum(total_infected_per_deme[1,:])
        maximum_infected_deme_2 = maximum(total_infected_per_deme[2,:])

        survived = total_infected[end] > 0
        survived_flag = survived ? 1 : 0
        extinction_time = !survived ? findfirst(total_infected .== 0) : NaN
        
        return (survived_flag, extinction_time, maximum_infected_deme_1, maximum_infected_deme_2)
    catch e
        println("error reading file: error is $e")
        return (0, NaN, NaN, NaN) # Assuming default values in case of error
    end
end

function analyze_migration_rates(MIGRATION_RATES, OUTPUT_DIRECTORY::String)
    survival_probabilities = fill(NaN, length(MIGRATION_RATES))

    # Your existing process_replicate function remains unchanged

    for (idx, migration_rate) in enumerate(MIGRATION_RATES)
        output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$idx")
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
        maximum_infected_deme_2 = [result[4] for result in results]

        valid_extinction_times = extinction_times[.!isnan.(extinction_times)]

        println("Maximum Extinction Time: ", maximum(valid_extinction_times))
        println("Average Extinction Time: ", mean(valid_extinction_times))
        println("number of times greater than 60: ", count(valid_extinction_times .> 60))
        println("number times greater than 70: ", count(valid_extinction_times .> 70))
        println("number times greater than 80: ", count(valid_extinction_times .> 80))
        println("Peak infected in deme 1: ", mean(maximum_infected_deme_1), " +- ", std(maximum_infected_deme_1))
        println("Peak infected in deme 2: ", mean(maximum_infected_deme_2), " +- ", std(maximum_infected_deme_2))

        survival_probabilities[idx] = sum(survived_results) / num_replicates
    end

    return survival_probabilities
end

survival_probabilities = analyze_migration_rates(MIGRATION_RATES, OUTPUT_DIRECTORY)

# Process and save trajectories for migration rate index 1
process_trajectories(1, 200)

# Convert to DataFrame for CSV
df = DataFrame(
    MigrationRate = MIGRATION_RATES,
    SurvivalProbability = survival_probabilities
)

# Save to CSV
CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "analysis_results.csv"), df)

# Process and save histograms for each migration rate, except the first and last
for idx in 2:length(MIGRATION_RATES)-1
    println(idx)
    flush(stdout)
    process_and_save_histograms(idx)
end

function process_special_case(migration_rate_idx)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$migration_rate_idx")
    replicate_files = glob("*.jld2", output_subdirectory)

    antigenic_variance_deme1 = []
    survival_status = []

    for file in replicate_files
        total_infected_per_deme, antigenic_variance_per_deme = read_data(file)

        # Find index of maximal infection for deme 1
        max_infected_idx_deme1 = argmax(total_infected_per_deme[1, :])

        # Record the antigenic variance at this point
        push!(antigenic_variance_deme1, antigenic_variance_per_deme[1, max_infected_idx_deme1])

        # Record survival status
        survived = sum(total_infected_per_deme[:, end]) > 0
        push!(survival_status, survived)
    end

    # Save the antigenic variance data along with survival status
    antigenic_variance_df = DataFrame(AntigenicVariance=antigenic_variance_deme1, Survival=survival_status)
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "special_case_antigenic_variance_deme1_migration_rate_idx_$migration_rate_idx.csv"), antigenic_variance_df)
end

# Call the function for migration_rate_idx = 1
process_special_case(1)

