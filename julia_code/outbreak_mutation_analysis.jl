using Glob
using Serialization
using Statistics
using Plots
using Printf
using Base.Threads
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

# const DIRECTORY_PATH = "simresults_newseed/"

const DIRECTORY_PATH ="C:/Users/Daniel/Desktop/simresults/"
const OUTPUT_PATH = "plotted_results_new/"

const MIGRATION_RATES = exp10.(LinRange(-7.0, -0.5, 9)) # Example migration rates to sweep over

function solve_S(R0::Float64; tol=1e-6, max_iter=1000)
    S = 0.5  # initial guess
    for _ in 1:max_iter
        new_S = exp(-R0 * (1-S))
        if abs(new_S - S) < tol
            break
        end
        S = new_S
    end
    return S
end

function calculate_secondary_infected_per_deme(simulation::Simulation, xstar::Float64)
    # Get the number of populations (demes)
    num_demes = length(simulation.trajectory[1].populations)
    num_time_points = length(simulation.duration_times)
    
    # Initialize a list of vectors to hold the total number of infected individuals at each time step for each deme
    secondary_infected_per_deme = zeros(num_demes, num_time_points);

    antigentic_region = abs.(simulation.trajectory[1].populations[1].xs) .> xstar

    # Iterate over each network state in the simulation's trajectory
    for (i, network) in enumerate(simulation.trajectory)
        # Iterate over each population in the current network state
        for (j, population) in enumerate(network.populations)
            # Increment the total infected count for the current time step and deme
            # by adding the integral of the viral density (viral_density multiplied by dx)
            secondary_infected_per_deme[j,i] += sum(population.viral_density[antigentic_region] .* population.dx)
        end
    end

    # Return the list containing the total number of infected individuals at each time step for each deme
    return secondary_infected_per_deme
end

function check_if_trajectory_survives(simulation::Simulation; cutoff=1e3)
    return calculate_total_infected(simulation)[end] > cutoff
end

function calculate_outbreak_probabilities(migration_rate, xstar; cutoff = 100, time_error_margin = 2.0)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    total_files = length(files)
    println(total_files)

    # Counters for each outcome
    count_deme1_first = 0
    count_deme2_first = 0
    count_both_same_time = 0
    count_survived = 0  # Counter for survived trajectories

    @threads for file in files
        simulation = open(deserialize, file)
        if check_if_trajectory_survives(simulation)
            count_survived += 1
            secondary_infected_per_deme = calculate_secondary_infected_per_deme(simulation, xstar)

            # Determine outbreak timings
            time_crosses_cutoff_deme1 = findfirst(secondary_infected_per_deme[1, :] .> cutoff)
            time_crosses_cutoff_deme2 = findfirst(secondary_infected_per_deme[2, :] .> cutoff)
            
            # Update counters based on outbreak timings
            if !isnothing(time_crosses_cutoff_deme1) && !isnothing(time_crosses_cutoff_deme2)
                if abs(time_crosses_cutoff_deme1 - time_crosses_cutoff_deme2) <= time_error_margin
                    count_both_same_time += 1
                elseif time_crosses_cutoff_deme1 < time_crosses_cutoff_deme2
                    count_deme1_first += 1
                else
                    count_deme2_first += 1
                end
            elseif !isnothing(time_crosses_cutoff_deme1)
                count_deme1_first += 1
            elseif !isnothing(time_crosses_cutoff_deme2)
                count_deme2_first += 1
            end
        end
    end
    
    # Calculate probabilities
    prob_deme1_first = count_deme1_first / count_survived
    prob_deme2_first = count_deme2_first / count_survived
    prob_both_same_time = count_both_same_time / count_survived

    error_deme1 = sqrt(prob_deme1_first * (1 - prob_deme1_first) / count_survived)
    error_deme2 = sqrt(prob_deme2_first * (1 - prob_deme2_first) / count_survived)
    error_both = sqrt(prob_both_same_time * (1 - prob_both_same_time) / count_survived)

    return prob_deme1_first, error_deme1, prob_deme2_first, error_deme2, prob_both_same_time, error_both

end

# Initialize arrays to store probabilities
prob_deme1_first_rates = zeros(length(MIGRATION_RATES))
prob_deme2_first_rates = zeros(length(MIGRATION_RATES))
prob_both_same_time_rates = zeros(length(MIGRATION_RATES))

error_deme1_rates = zeros(length(MIGRATION_RATES))
error_deme2_rates = zeros(length(MIGRATION_RATES))
error_both_rates = zeros(length(MIGRATION_RATES))

# Loop through migration rates
for (i, migration_rate) in enumerate(MIGRATION_RATES)
    println(i)
    # Calculate probabilities for each migration rate
    prob_deme1_first, error_deme1, prob_deme2_first, error_deme2, prob_both_same_time, error_both = calculate_outbreak_probabilities(migration_rate, 2.7; cutoff = 1000, time_error_margin = 2.0)

    # Store the calculated probabilities and errors
    prob_deme1_first_rates[i] = prob_deme1_first
    prob_deme2_first_rates[i] = prob_deme2_first
    prob_both_same_time_rates[i] = prob_both_same_time

    error_deme1_rates[i] = error_deme1
    error_deme2_rates[i] = error_deme2
    error_both_rates[i] = error_both
end

# Create a plot
p = plot(xscale=:log10, xlabel="Migration Rate", ylabel="Probability", title="Probability of Outbreak Occurrence by Migration Rate", legend=:topleft)

# Ribbon plots
plot!(p, MIGRATION_RATES, prob_deme1_first_rates, ribbon=error_deme1_rates, label="Deme 1 First", marker=:circle)
plot!(p, MIGRATION_RATES, prob_deme2_first_rates, ribbon=error_deme2_rates, label="Deme 2 First", marker=:square)
plot!(p, MIGRATION_RATES, prob_both_same_time_rates, ribbon=error_both_rates, label="Both at Same Time", marker=:triangle)

# Display the plot
display(p)

# Optionally, save the plot to a file
savefig(p, "$(OUTPUT_PATH)/spatial_distribution_secondary_outbreak.png")

# Create a plot
p = plot(xscale=:log10, xlabel="Migration Rate", ylabel="Probability", title="Probability of Outbreak Occurrence by Migration Rate", legend=:topright)

# Plot the areas between the curves
plot!(p, MIGRATION_RATES, prob_both_same_time_rates, fillrange=0, fillalpha=0.5, label="Both at Same Time", color=:green)
plot!(p, MIGRATION_RATES, prob_deme2_first_rates .+ prob_both_same_time_rates, fillrange=prob_both_same_time_rates, fillalpha=0.5, label="Deme 2 First", color=:red)
plot!(p, MIGRATION_RATES, prob_deme1_first_rates .+ prob_both_same_time_rates .+ prob_deme2_first_rates, fillrange=prob_both_same_time_rates .+ prob_deme2_first_rates, fillalpha=0.5, label="Deme 1 First", color=:blue)

# Display the plot
display(p)

# Optionally, save the plot to a file
savefig(p, "$(OUTPUT_PATH)/cumulative_outbreak_probabilities_migration_rate.png")
