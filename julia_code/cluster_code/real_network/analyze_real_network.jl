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

function process_and_save_histograms(outbreak_deme_idx; cutoff = 100)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "outbreak_deme_idx_$outbreak_deme_idx")
    replicate_files = glob("*.jld2", output_subdirectory)

    antigenic_variance_deme1 = []
    antigenic_variance_deme2 = []
    peak_time_difference = []
    variance_difference = []

    for file in replicate_files
        total_infected_per_deme, antigenic_variance_per_deme = read_data(file)
        
        # Find indices of maximal infection for each deme because dt = 1/20 and thin by = 20 index corrresponds to "real" time!
        max_infected_idx_deme1 = argmax(total_infected_per_deme[1, :])
        min_infected_idx_deme1 = argmin(total_infected_per_deme[1, :])
        
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
                println("LATE")
                flush(stdout)
            end
     
            if total_infected_per_deme[2, max_infected_idx_deme2] > cutoff && first_cross_idx_deme2 <= min_infected_idx_deme1
                push!(antigenic_variance_deme2, antigenic_variance_per_deme[2, max_infected_idx_deme2])
                push!(peak_time_difference, max_infected_idx_deme2 - max_infected_idx_deme1)
                push!(variance_difference, antigenic_variance_per_deme[2, max_infected_idx_deme2] - antigenic_variance_per_deme[1, max_infected_idx_deme1])
            end
        end
    end

    # Save the antigenic variance data
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "antigenic_variance_deme1_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(AntigenicVariance=antigenic_variance_deme1))
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "antigenic_variance_deme2_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(AntigenicVariance=antigenic_variance_deme2))
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "peak_time_difference_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(PeakTimeDifference=peak_time_difference))
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "variance_difference_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(VarianceDifference=variance_difference))
end

function process_trajectories(migration_rate_idx, num_replicates)
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


# Initialize survival_probabilities as an array with NaN
survival_probabilities = fill(NaN, length(MIGRATION_RATES))

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




