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
        
        # Initialize survival count
        survival_count = 0

        # Iterate over each replicate
        for replicate_data in values(all_replicates_data)
            # Check if the data structure contains total infected per deme
            if haskey(replicate_data, "total_infected")
                # Sum the total infected across all demes for the last time point
                total_infected_sum = sum([deme_data[end] for deme_data in replicate_data["total_infected"]])

                # Update survival count if total infected sum is greater than 0
                if total_infected_sum > 0
                    survival_count += 1
                end
            end
        end

        # Calculate probability of survival
        probability_of_survival = survival_count / num_replicates
        return probability_of_survival, num_replicates
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

        # Initialize a plot with three subplots
        p1 = plot(;yscale = :log10, ylims = (1, Inf)) # Total infected vs time
        p2 = plot(;yscale = :log10, ylims = (1, Inf)) # Infected in Deme 1 vs time
        p3 = plot(;yscale = :log10, ylims = (1, Inf)) # Infected in all demes but Deme 1 vs time

        for (replicate_idx, data) in pairs(all_replicates_data)
            total_infected = data["total_infected"]
            times = data["times"]

            # Sum over all demes for total infected
            total_infected_sum = sum(hcat(total_infected...), dims=2)
            reg = total_infected_sum .> 0
            plot!(p1, times[vec(reg)], vec(total_infected_sum[reg]), label=false)

            # Plot for Deme 1
            reg = total_infected[1] .> 0
            plot!(p2, times[reg], total_infected[1][reg], label=false)

            # Sum over all demes except Deme 1
            total_infected_except_deme1 = sum(hcat(total_infected[2:end]...), dims=2)
            reg = total_infected_except_deme1 .> 0
            plot!(p3, times[vec(reg)], vec(total_infected_except_deme1[reg]), label=false)
        end

        # Customize the subplots
        xlabel!(p1, "Time"); ylabel!(p1, "Total Infected")
        xlabel!(p2, "Time"); ylabel!(p2, "Infected in Deme 1")
        xlabel!(p3, "Time"); ylabel!(p3, "Infected in All Demes Except Deme 1")

        # Combine subplots into one figure with increased size
        combined_plot = plot(p1, p2, p3, layout = (3, 1), title = ["Total Infected" "Infected in Deme 1" "Infected in All Demes Except Deme 1"], size = (600, 900))

        # Save the plot to file
        savefig(joinpath(PLOT_RESULTS_DIRECTORY, "Infected_Trajectories_N$(network_size)_M$(migration_rate_idx).png"))
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
