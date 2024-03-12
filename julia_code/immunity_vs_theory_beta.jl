using Glob
using Serialization
using Statistics
using Plots
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DIRECTORY_PATH ="C:/Users/Daniel/Desktop/simresults_varying_beta/"
const OUTPUT_PATH = "./plotted_results_beta_sweep/"

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

function total_immunity_at_time(simulation::CoevolutionNetworkBase.Simulation, time_target::Float64)
    # Find index closest to target time
    idx_target = argmin(abs.(simulation.duration_times .- time_target))

    return sum(simulation.trajectory[idx_target].populations[1].immune_density .* simulation.trajectory[idx_target].populations[1].dx)
end

function empirical_total_immunity_for_beta(beta_value, time_target::Float64, alpha, gamma, M)
    files = glob("simulation_results_beta_$(beta_value)_replicate_*.jld2", DIRECTORY_PATH)
    immunities = zeros(Float64, length(files))

    Threads.@threads for idx in 1:length(files)
        file = files[idx]
        simulation = open(deserialize, file)
        
        
        immunities[idx] = total_immunity_at_time(simulation, time_target)
        
    end

    reg = immunities .> 0.05
    avg_immunity = mean(immunities[reg])
    std_error_immunity = std(immunities[reg]) / sqrt(sum(reg))

    return avg_immunity, std_error_immunity
end

function main()
    betas = LinRange(1.1, 6, 10)
    time_target = 24.0

    alpha = 1.0  # sample value, modify as needed
    gamma = 0.0  # sample value, modify as needed
    M = 15       # sample value, modify as needed

    empirical_immunities = zeros(Float64, length(betas))
    empirical_errors = zeros(Float64, length(betas))
    theoretical_immunities = zeros(Float64, length(betas))

    for idx in 1:length(betas)
        beta_value = betas[idx]
        println(beta_value)
        
        avg_immunity, std_error = empirical_total_immunity_for_beta(beta_value, time_target, alpha, gamma, M)
        empirical_immunities[idx] = avg_immunity
        empirical_errors[idx] = std_error

        R0 = beta_value / (alpha + gamma)
        S = solve_S(R0)
        theoretical_immunities[idx] = 1 - (S)^(1/M)
    end

    # Plot the empirical total immunity with error bars and the theoretical values
    plot(betas, empirical_immunities, ribbon=empirical_errors, label="Empirical", legend=:topright, 
         xlabel="Beta", ylabel="Total Immunity at time $time_target", alpha=0.5)
    plot!(betas, theoretical_immunities, label="Theoretical", color="red")

    # Saving the plot
    savefig(joinpath(OUTPUT_PATH, "total_immunity_vs_beta.png"))
end

main()
