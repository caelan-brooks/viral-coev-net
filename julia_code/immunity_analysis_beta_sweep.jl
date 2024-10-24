using Glob
using Serialization
using Statistics
using Plots
using Printf
using Base.Threads
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DIRECTORY_PATH ="C:/Users/Daniel/Desktop/simresults_varying_beta/"
const OUTPUT_PATH = "plotted_results_beta_sweep/"

function calculate_total_immunity_per_deme(simulation::CoevolutionNetworkBase.Simulation)
    num_demes = length(simulation.trajectory[1].populations)
    total_immunity_per_deme = [zeros(Float64, length(simulation.trajectory)) for _ in 1:num_demes]

    for (i, network) in enumerate(simulation.trajectory)
        for (j, population) in enumerate(network.populations)
            total_immunity_per_deme[j][i] += sum(population.immune_density .* population.dx)
        end
    end
    return total_immunity_per_deme
end

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

function plot_total_immunity_per_deme_for_beta(beta_value)
    files = glob("simulation_results_beta_$(beta_value)_replicate_*.jld2", DIRECTORY_PATH)
    formatted_beta_value = @sprintf("%.2e", beta_value)

    println("Number of files for beta = $formatted_beta_value: ", length(files))

    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    # Fetch parameters from the first simulation (assuming they are consistent across files)
    sample_simulation = open(deserialize, files[1])
    alpha = sample_simulation.trajectory[1].populations[1].alpha
    gamma = sample_simulation.trajectory[1].populations[1].gamma
    M = sample_simulation.trajectory[1].populations[1].M

    R0 = beta_value / (alpha + gamma)
    S = solve_S(R0)
    theoretical_immunity = 1 - (S)^(1/M)

    # Data structure to store results from threads
    results = Vector{Tuple{Vector{Float64}, Vector{Vector{Float64}}}}(undef, length(files))
    
    # Multi-threaded data calculation
    Threads.@threads for idx in 1:length(files)
        file = files[idx]
        simulation = open(deserialize, file)
        xdata = simulation.duration_times
        ydata_list = calculate_total_immunity_per_deme(simulation)
        
        results[idx] = (xdata, ydata_list)
    end

    # Single-threaded plotting
    plots_immunity = [plot(title="Total Immunity (Beta: $formatted_beta_value, Deme: $i)", 
                           xlabel="Time", 
                           ylabel="Total Immunity in Deme $i", 
                           legend=false) for i in 1:2]
    
    for (xdata, ydata_list) in results
        for (idx, ydata) in enumerate(ydata_list)
            plot!(plots_immunity[idx], xdata, ydata, label=false, legend=false)
        end
    end

    # Plot theoretical immunity on top of each curve
    for p in plots_immunity
        hline!(p, [theoretical_immunity], linestyle=:dash, linewidth=2, color=:black)
    end

    for (idx, p) in enumerate(plots_immunity)
        savefig(p, "$(OUTPUT_PATH)/total_immunity_deme_$(idx)_beta_$(formatted_beta_value).png")
    end
end

BETA_VALUES = LinRange(1.1, 6, 10)  # Example beta values to sweep over
for beta in BETA_VALUES
    plot_total_immunity_per_deme_for_beta(beta)
end
