using Glob
using Serialization
using Plots
using Printf
using Base.Threads
using Statistics

include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DIRECTORY_PATH = "C:/Users/Daniel/Desktop/simresults_oneinfected/"
const OUTPUT_PATH = "plotted_results_oneinfected"

function calculate_total_infected_per_deme(simulation::Simulation)
    # Get the number of populations (demes)
    num_demes = length(simulation.trajectory[1].populations)
    
    # Initialize a list of vectors to hold the total number of infected individuals at each time step for each deme
    total_infected_per_deme = [zeros(Float64, length(simulation.trajectory)) for _ in 1:num_demes]

    # Iterate over each network state in the simulation's trajectory
    for (i, network) in enumerate(simulation.trajectory)
        # Iterate over each population in the current network state
        for (j, population) in enumerate(network.populations)
            # Increment the total infected count for the current time step and deme
            # by adding the integral of the viral density (viral_density multiplied by dx)
            total_infected_per_deme[j][i] += sum(population.viral_density .* population.dx)
        end
    end

    # Return the list containing the total number of infected individuals at each time step for each deme
    return total_infected_per_deme
end

function calculate_prob_above_cutoff(migration_rate, cutoff)
    # Get the list of files for the given migration rate
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    # Counter for the number of replicates where deme 2 exceeds the cutoff
    count_exceeding_cutoff = 0
    
    for file in files
        simulation = open(deserialize, file)
        
        # Get the total infected for each deme over time
        ydata_list = calculate_total_infected_per_deme(simulation)

        # Only consider deme 2 data
        deme_2_data = ydata_list[2]

        # Check if the maximum population in deme 2 is above the cutoff
        if maximum(deme_2_data) > cutoff
            count_exceeding_cutoff += 1
        end
    end
    
    # Calculate the probability
    probability = count_exceeding_cutoff / length(files)
    return probability
end

function main_task_1()
    migration_rates = vcat([0], exp10.(LinRange(-6, 0.5, 9))) # Example migration rates to sweep over
    cutoff = 10^5

    # Store probabilities for each migration rate
    probabilities = Float64[]

    for migration_rate in migration_rates
        println(migration_rate)
        # Skip migration rate of 0 as specified
        if migration_rate == 0
            push!(probabilities, NaN)
            continue
        end
        prob = calculate_prob_above_cutoff(migration_rate, cutoff)
        push!(probabilities, prob)
    end
    
    # Plot the results
    p = plot(migration_rates[2:end], probabilities[2:end], xscale=:log10, xlabel="Migration Rate", ylabel="Probability",
         title="Probability of Deme 2 Population Exceeding $cutoff", legend=false, marker=:circle)

    # Check if the "plotted_results" folder exists in the script's directory, if not create it
    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    # Save the plot
    savefig(p, "$(OUTPUT_PATH)/probability_exceeding_cutoff.png")
end

print("Starting task 1...\n")
# main_task_1()

function calculate_time_to_max(migration_rate, cutoff)
    # Get the list of files for the given migration rate
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    # List to store the times it took to reach the maximum population for each replicate
    times_to_max = Float64[]

    for file in files
        simulation = open(deserialize, file)
        
        # Get the total infected for each deme over time
        ydata_list = calculate_total_infected_per_deme(simulation)

        # Only consider deme 2 data
        deme_2_data = ydata_list[2]

        # Check if the maximum population in deme 2 is above the cutoff
        if maximum(deme_2_data) > cutoff
            # Find the time when it reaches the maximum
            time_idx = findfirst(isequal(maximum(deme_2_data)), deme_2_data)
            push!(times_to_max, simulation.duration_times[time_idx])
        end
    end

    # Return average time to max if there are any valid times, otherwise return NaN
    return isempty(times_to_max) ? NaN : mean(times_to_max)
end

function main_task_2()
    migration_rates = vcat([0], exp10.(LinRange(-6, 0.5, 9))) # Example migration rates to sweep over
    cutoff = 10^5

    # Store average times for each migration rate
    avg_times = Float64[]

    for migration_rate in migration_rates
        println(migration_rate)
        avg_time = calculate_time_to_max(migration_rate, cutoff)
        push!(avg_times, avg_time)
    end
    
    # Plot the results
    p = plot(migration_rates[2:end], avg_times[2:end], xscale=:log10, xlabel="Migration Rate", ylabel="Average Time to Max Population",
         title="Average Time for Deme 2 Population to Reach Max Above $cutoff", legend=false, marker=:circle)

    # Check if the "plotted_results" folder exists in the script's directory, if not create it
    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    # Save the plot
    savefig(p, "$(OUTPUT_PATH)/avg_time_to_max_population.png")
end

print("Starting task 2...\n")
# main_task_2()

function calculate_time_difference_to_max(migration_rate, cutoff)
    # Get the list of files for the given migration rate
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    # List to store the time differences for each replicate
    time_diffs = Float64[]

    Threads.@threads for file in files
        simulation = open(deserialize, file)
        
        # Get the total infected for each deme over time
        ydata_list = calculate_total_infected_per_deme(simulation)

        # Consider both deme 1 and 2 data
        deme_1_data = ydata_list[1]
        deme_2_data = ydata_list[2]

        # Check if the maximum population in deme 2 is above the cutoff
        if maximum(deme_2_data) > cutoff
            # Find the time when each deme reaches the maximum
            time_idx_deme_1 = findfirst(isequal(maximum(deme_1_data)), deme_1_data)
            time_idx_deme_2 = findfirst(isequal(maximum(deme_2_data)), deme_2_data)
            push!(time_diffs, simulation.duration_times[time_idx_deme_2] - simulation.duration_times[time_idx_deme_1])
        end
    end

    # Return average time difference if there are any valid times, otherwise return NaN
    return isempty(time_diffs) ? NaN : mean(time_diffs)
end

using LinearAlgebra

function power_law_fit(x, y)
    # Ensure data points are non-zero and filter out any that aren't
    valid_points = (x .> 0) .& (y .> 0)
    x = x[valid_points]
    y = y[valid_points]

    # Convert data to log-log scale
    logx = log.(x)
    logy = log.(y)

    # Compute the linear regression coefficients
    A = [ones(length(logx)) logx]
    coef = A \ logy

    # The slope in log-log scale will give us the power-law exponent
    b = coef[2]
    a = exp(coef[1])

    return a, b
end

function main_task_3()
    migration_rates = vcat([0], exp10.(LinRange(-6, 0.5, 9))) # Example migration rates to sweep over
    cutoff = 10^5

    # Store average time differences for each migration rate
    avg_time_diffs = Float64[]

    for migration_rate in migration_rates
        println(migration_rate)
        avg_time_diff = calculate_time_difference_to_max(migration_rate, cutoff)
        push!(avg_time_diffs, avg_time_diff)
    end
    
    # Plot the results
    reg = (migration_rates .> 0) .& (avg_time_diffs .> 0)
    p = plot(migration_rates[reg], avg_time_diffs[reg], xscale=:log10, yscale=:log10, xlabel="Migration Rate", ylabel="Average Time Difference",
         title="Average Time Difference for Deme 2 and Deme 1 to Reach Max Above $cutoff", legend=true, marker=:circle)


    # Fit the first four points to a power law
    x_vals = migration_rates[reg][1:6]
    y_vals = avg_time_diffs[reg][1:6]
    a, b = power_law_fit(x_vals, y_vals)

    # Generate fitted y-values using the power-law relationship
    fitted_y_vals = a .* x_vals .^ b
    plot!(p, x_vals, fitted_y_vals, linestyle=:dash, color=:red, label=@sprintf("y = %.3f x^%.3f", a, b))

    # Check if the "plotted_results" folder exists in the script's directory, if not create it
    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    # Save the plot
    savefig(p, "$(OUTPUT_PATH)/avg_time_difference_to_max_population.png")
end

print("Starting task 3...\n")
main_task_3()