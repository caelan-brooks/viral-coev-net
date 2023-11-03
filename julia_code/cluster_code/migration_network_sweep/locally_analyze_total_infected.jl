using Serialization
using Glob
using Plots
using Statistics

# Constants
const LOCAL_RESULTS_DIRECTORY = "./migration_network_sweep_results/total_infected"
const PLOT_RESULTS_DIRECTORY = "./plotted_results"
const MIGRATION_RATES = exp10.(LinRange(-6, 0.5, 10))
const NETWORK_SIZES = 2:6

# Ensure the plot results directory exists
isdir(PLOT_RESULTS_DIRECTORY) || mkpath(PLOT_RESULTS_DIRECTORY)

# Function to calculate the probability of survival and the number of replicates
function calculate_probability_of_survival_and_replicates(network_size, migration_rate_idx)
    file_path = joinpath(LOCAL_RESULTS_DIRECTORY, "network_size_$(network_size)_migration_rate_idx_$(migration_rate_idx).jld2")
    if isfile(file_path)
        all_replicates_data = open(deserialize, file_path)
        num_replicates = length(all_replicates_data)
        survival_count = sum([!isempty(v["total_infected"]) && last(v["total_infected"]) > 0 for v in values(all_replicates_data)])
        return (survival_count / num_replicates), num_replicates
    else
        return NaN, 0 # File does not exist or no replicates
    end
end

# Function to calculate error for plotting
function calculate_error(probability, num_replicates)
    # Standard error of the mean for the binomial proportion
    return sqrt((probability * (1 - probability)) / num_replicates)
end

# Plot survival probabilities
function plot_survival_probabilities()
    plot() # Initialize plot

    for network_size in NETWORK_SIZES
        probabilities = Float64[]
        errors = Float64[]
        
        for idx in 1:length(MIGRATION_RATES)
            prob, num_reps = calculate_probability_of_survival_and_replicates(network_size, idx)
            push!(probabilities, prob)
            if num_reps > 0
                push!(errors, calculate_error(prob, num_reps))
            else
                push!(errors, NaN) # If no replicates, error is undefined
            end
        end

        # Add the probability line and ribbon
        plot!(MIGRATION_RATES, probabilities, ribbon=errors, label="Network Size $network_size", xaxis=:log)
    end

    # Customize and save the plot
    xlabel!("Migration Rate")
    ylabel!("Probability of Survival")
    title!("Survival Probability vs. Migration Rate")
    savefig(joinpath(PLOT_RESULTS_DIRECTORY, "Survival_Probability.png"))
end

# Function to plot total infected trajectories for a given network size and migration rate
function plot_total_infected_trajectories(network_size, migration_rate_idx)
    file_path = joinpath(LOCAL_RESULTS_DIRECTORY, "network_size_$(network_size)_migration_rate_idx_$(migration_rate_idx).jld2")
    
    if isfile(file_path)
        all_replicates_data = open(deserialize, file_path)

        p = plot(;yscale = :log10) # Initialize a new plot

        for (replicate_idx, data) in pairs(all_replicates_data)
            total_infected = data["total_infected"]
            times = data["times"]
            reg = total_infected .> 0
            plot!(p, times[reg], total_infected[reg], label=false) # Add each trajectory to the plot without a legend entry
        end

        # Customize the plot
        xlabel!(p, "Time")
        ylabel!(p, "Total Infected")
        title!(p, "Infected Trajectories\nNetwork Size: $network_size, Migration Rate Index: $migration_rate_idx")

        # Save the plot to file
        savefig(p, joinpath(PLOT_RESULTS_DIRECTORY, "Infected_Trajectories_N$(network_size)_M$(migration_rate_idx).png"))
    else
        println("File not found: $file_path")
    end
end

# Function to plot total infected trajectories for all network sizes and migration rates
function plot_all_infected_trajectories()
    for network_size in NETWORK_SIZES
        for migration_rate_idx in 1:length(MIGRATION_RATES)
            plot_total_infected_trajectories(network_size, migration_rate_idx)
        end
    end
end

# Modify the main function to include the new plotting function
function main()
    plot_survival_probabilities()
    plot_all_infected_trajectories()
end

main()
