using Serialization
using CSV
using DataFrames
using Glob

const OUTPUT_DIRECTORY = "/pool001/dswartz/three_deme_one_edge"
const CSV_OUTPUT_DIRECTORY = "/pool001/dswartz/three_deme_one_edge/csv_outputs"  # Directory for CSV outputs
const MIGRATION_RATES = [0; exp10.(LinRange(-7, -0.5, 10))]

# Create CSV output directory if it doesn't exist
isdir(CSV_OUTPUT_DIRECTORY) || mkdir(CSV_OUTPUT_DIRECTORY)

# Define read_data as a global function
function read_data(file)
    open(deserialize, file)
end

function process_and_save_histograms(migration_rate_idx)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$migration_rate_idx")
    replicate_files = glob("*.jld2", output_subdirectory)

    antigenic_variance_deme1 = []
    antigenic_variance_deme2 = []

    for file in replicate_files
        total_infected_per_deme, antigenic_variance_per_deme = read_data(file)
        
        # Find indices of maximal infection for each deme
        max_infected_idx_deme1 = argmax(total_infected_per_deme[1, :])
        max_infected_idx_deme2 = argmax(total_infected_per_deme[2, :])

        # Record the antigenic variance at these points
        push!(antigenic_variance_deme1, antigenic_variance_per_deme[1, max_infected_idx_deme1])
        push!(antigenic_variance_deme2, antigenic_variance_per_deme[2, max_infected_idx_deme2])
    end

    # Save the antigenic variance data
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "antigenic_variance_deme1_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(AntigenicVariance=antigenic_variance_deme1))
    CSV.write(joinpath(CSV_OUTPUT_DIRECTORY, "antigenic_variance_deme2_migration_rate_idx_$migration_rate_idx.csv"), DataFrame(AntigenicVariance=antigenic_variance_deme2))
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
for (idx, migration_rate) in enumerate(MIGRATION_RATES)
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$idx")
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
    # process_and_save_histograms(idx)
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
# process_special_case(1)

