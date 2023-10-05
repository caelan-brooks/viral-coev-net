using Glob
using Serialization
using Statistics
using Plots
using Printf
using Base.Threads
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DIRECTORY_PATH ="C:/Users/Daniel/Desktop/simresults_oneinfected_larger_beta/"
const OUTPUT_PATH = "plotted_results_oneinfected_larger_beta/"

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

function plot_total_immunity_per_deme(migration_rate, migration_rate_index)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    formatted_migration_rate = @sprintf("%.2e", migration_rate)

    println("Number of files: ", length(files))

    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    # Fetch parameters from the first simulation (assuming they are consistent across files)
    sample_simulation = open(deserialize, files[1])
    beta = sample_simulation.trajectory[1].populations[1].beta
    alpha = sample_simulation.trajectory[1].populations[1].alpha
    gamma = sample_simulation.trajectory[1].populations[1].gamma
    M = sample_simulation.trajectory[1].populations[1].M

    R0 = beta / (alpha + gamma)
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
    plots_immunity = [plot(title="Total Immunity (Migration Rate: $formatted_migration_rate, Deme: $i)", 
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
        savefig(p, "$(OUTPUT_PATH)/total_immunity_deme_$(idx)_idx_$(migration_rate_index)_rate_$(formatted_migration_rate).png")
    end
end


migration_rates = vcat([0], exp10.(LinRange(-6, 0.5, 9)))
for (index, rate) in enumerate(migration_rates)
    plot_total_immunity_per_deme(rate, index)
end

