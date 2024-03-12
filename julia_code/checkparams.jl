using Glob
using Serialization
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const BASE_DIRECTORY = "C:/Users/Daniel/Desktop/simresults_oneinfected"

# Possible values for "SOMETHING"
possible_values = ["", "_large_r", "_larger_beta","_smaller_mutation"]  # Replace with actual values

function print_population_data_without_densities(directory)
    # Get all files in the directory
    files = glob("*.jld2", directory)
    
    # If there are no files in the directory, return
    if isempty(files)
        println("No files found in $directory")
        return
    end

    # Pick the first file (it doesn't matter which one)
    file_path = files[1]
    
    # Deserialize the simulation from the file
    simulation = open(deserialize, file_path)
    
    # Extract a single population object (assuming the population object is stored directly)
    # Modify as necessary based on the structure of your data
    population = simulation.trajectory[1].populations[1]  # Modify this line based on how your data is structured
    
    # Print the population's data except for `viral_density`, `immune_density`, and the other specified fields
    for field in fieldnames(typeof(population))
        if field ∉ (:viral_density, :immune_density, :cross_reactive, :fitness, :susceptibility, :temporary_data, :xs)
            println("$field: ", getfield(population, field))
        end
    end

end

# Iterate through each possible value and execute the function
for value in possible_values
    directory = BASE_DIRECTORY * value
    println("\nProcessing directory: $directory\n")
    print_population_data_without_densities(directory)
end
