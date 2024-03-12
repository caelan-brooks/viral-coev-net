using Serialization
using Glob
using Printf
using Statistics
using Plots
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DATA_DIRECTORY = "C:/Users/Daniel/Desktop/simresults_varying_beta"
const OUTPUT_PATH = "plotted_results_beta_sweep/"
const BETA_VALUES = LinRange(1.1, 6, 10)  # Example beta values to sweep over

"""
    calculate_probability_of_survival(beta_value, cutoff)

Calculate the probability of survival for a given beta value.
"""
function calculate_probability_of_survival(beta_value, cutoff)
    files = glob("simulation_results_beta_$(beta_value)_replicate_*.jld2", DATA_DIRECTORY)
    println(length(files))
    
    survival_counts = 0
    total_files = length(files)
    survival_outcomes = zeros(Int,length(files))

    Threads.@threads for idx in 1:length(files)
        file = files[idx]
        try
            simulation = open(deserialize, file)
            data_slice = calculate_total_infected(simulation)[Int(round(end * 0.95)) : end]

            if mean(data_slice) > cutoff
                survival_outcomes[idx] += 1
            end
        catch e
            println("Error processing file $file: $e")
        end
    end
    
    probability_of_survival = sum(survival_outcomes) / total_files
    standard_error = sqrt((probability_of_survival * (1 - probability_of_survival)) / total_files)
    
    return probability_of_survival, standard_error
end

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

function plot_total_infected_per_deme_for_beta(beta_value)
    files = glob("simulation_results_beta_$(beta_value)_replicate_*.jld2", DATA_DIRECTORY)
    formatted_beta_value = @sprintf("%.2e", beta_value)
    total_files = length(files)

    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    # Data structure to store results from threads
    processed_data = Vector{Tuple{Vector{Float64}, Vector{Vector{Float64}}, Vector{Float64}}}(undef, total_files)

    # Multi-threaded data extraction and processing
    Threads.@threads for idx in 1:total_files
        file = files[idx]
        simulation = open(deserialize, file)
        xdata = simulation.duration_times
        ydata_list = calculate_total_infected_per_deme(simulation)
        ydata_end_values = [ydata[end] for ydata in ydata_list]

        processed_data[idx] = (xdata, ydata_list, ydata_end_values)
    end

    # Initialize plots
    plots = [plot(title="Trajectories (Beta: $formatted_beta_value, Deme: $i)", xlabel="Time", ylabel="Total Infected Number in Deme $i", legend=false, ylims=(1,3*10^6)) for i in 1:2]

    # Single-threaded plotting
    for (xdata, ydata_list, _) in processed_data
        for (idx, ydata) in enumerate(ydata_list)
            reg = ydata .> 0.0
            if any(reg)
                plot!(plots[idx], xdata[reg], ydata[reg], label=false, legend=false, yscale=:log10)
            else
                plot!(plots[idx], xdata[reg], ydata[reg], label=false, legend=false)
            end
        end
    end

    for (idx, p) in enumerate(plots)
        savefig(p, "$(OUTPUT_PATH)/trajectory_deme_$(idx)_beta_$(formatted_beta_value).png")
    end
end

function main()

    # Calculate probability of survival
    cutoff = 1000
    results = Dict{Float64, Tuple{Float64, Float64}}()
    betas = Float64[]
    survival_probs = Float64[]
    errors = Float64[]

    for beta_value in BETA_VALUES
        prob_survival, std_error = calculate_probability_of_survival(beta_value, cutoff)
        results[beta_value] = (prob_survival, std_error)
        
        # Store values for plotting
        push!(betas, beta_value)
        push!(survival_probs, prob_survival)
        push!(errors, std_error)

        println("Beta: $beta_value, Probability of Survival: $prob_survival ± $std_error")
    end

    # Plot survival probability with error bars
    p = plot(betas, survival_probs, ribbon=errors, fillalpha=0.3,
             xlabel="Beta", ylabel="Survival Probability",
             title="Survival Probability as a function of Beta",
             legend=false, linewidth=2, linecolor=:blue, yscale=:log10)
    scatter!(p, betas, survival_probs, color=:red, markersize=5)
    plot!(p, size=(800, 600)) # to set the size of the plot
   
    savefig(p, joinpath(OUTPUT_PATH, "survival_probability_vs_beta.png"))

    # Plot the trajectories
    for beta_value in BETA_VALUES
        plot_total_infected_per_deme_for_beta(beta_value)
    end
end

main()

