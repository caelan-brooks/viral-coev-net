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

function calculate_probabilities(beta_value, cutoff, drift_time_cutoff)
    files = glob("simulation_results_beta_$(beta_value)_replicate_*.jld2", DATA_DIRECTORY)

    drift_extinct_outcomes = zeros(Int, length(files))
    immunity_extinct_outcomes = zeros(Int, length(files))

    Threads.@threads for idx in 1:length(files)
        file = files[idx]
        simulation = open(deserialize, file)
        total_infected_end = calculate_total_infected(simulation)[end]
        total_infected_drift_time = calculate_total_infected(simulation)[Int(round(length(simulation.duration_times) * drift_time_cutoff / simulation.duration_times[end]))]

        # Check drift extinction: the trajectory goes extinct before drift_time_cutoff
        if total_infected_drift_time < cutoff
            drift_extinct_outcomes[idx] += 1
        elseif total_infected_end < cutoff # Check immunity extinction: survives drift but goes extinct by the end
            immunity_extinct_outcomes[idx] += 1
        end
    end

    prob_drift_extinct = sum(drift_extinct_outcomes) / length(files)
    prob_immunity_extinct = sum(immunity_extinct_outcomes) / length(files)
    prob_survival = 1 - (prob_drift_extinct + prob_immunity_extinct)

    # Compute standard errors
    se_drift_extinct = sqrt(prob_drift_extinct * (1 - prob_drift_extinct) / length(files))
    se_immunity_extinct = sqrt(prob_immunity_extinct * (1 - prob_immunity_extinct) / length(files))
    se_survival = sqrt(prob_survival * (1 - prob_survival) / length(files))

    return prob_drift_extinct, prob_immunity_extinct, prob_survival, se_drift_extinct, se_immunity_extinct, se_survival
end


function is_extinct_before_time(simulation, time_cutoff)
    idx_cutoff = argmin(abs.(simulation.duration_times .- time_cutoff))
    return sum(simulation.trajectory[idx_cutoff].populations[1].viral_density) == 0
end

function main()

    # Calculate probabilities
    cutoff = 10
    drift_time_cutoff = 5.0  # Example value, adjust as needed

    files = glob("simulation_results_beta_$(BETA_VALUES[1])_replicate_*.jld2", DATA_DIRECTORY)
    file = files[1]
    sample_simulation = open(deserialize,file)
    alpha = sample_simulation.trajectory[1].populations[1].alpha # Fetch from a sample simulation

    betas = []
    drift_extinct_probs = []
    immunity_extinct_probs = []
    survival_probs = []
    se_drift_extinct_list = []
    se_immunity_extinct_list = []
    se_survival_list = []
    
    for beta_value in BETA_VALUES
        println(beta_value)
        prob_drift_extinct, prob_immunity_extinct, prob_survival, se_drift_extinct, se_immunity_extinct, se_survival = calculate_probabilities(beta_value, cutoff, drift_time_cutoff)

        # Store values for plotting
        push!(betas, beta_value)
        push!(drift_extinct_probs, prob_drift_extinct)
        push!(immunity_extinct_probs, prob_immunity_extinct)
        push!(survival_probs, prob_survival)

        push!(se_drift_extinct_list, se_drift_extinct)
        push!(se_immunity_extinct_list, se_immunity_extinct)
        push!(se_survival_list, se_survival)

        println("Beta: $beta_value, 
                Drift Extinct Probability: $prob_drift_extinct, 
                Immunity Extinct Probability: $prob_immunity_extinct,
                Survival Probability: $prob_survival")
    end
    
    theoretical_drift_extinct_probs = alpha ./ betas # Theoretical drift extinction probability for each beta
    # Plotting
    p = plot(betas, drift_extinct_probs, ribbon=se_drift_extinct_list, fillalpha=0.3, label="Empirical Drift Extinction", legend=:topright)
    # plot!(p, betas, theoretical_drift_extinct_probs, line=:dash, color=:blue, label="Theoretical Drift Extinction")
    plot!(p, betas, immunity_extinct_probs, ribbon=se_immunity_extinct_list, fillalpha=0.3, label="Immunity Extinction")
    plot!(p, betas, survival_probs, ribbon=se_survival_list, fillalpha=0.3, label="Survival")
    
    # Format the plot
    plot!(p, xlabel="Beta", ylabel="Probability", title="Outcomes as a function of Beta", size=(800, 600))

    # Save the plot
    savefig(p, joinpath(OUTPUT_PATH, "outcomes_vs_beta.png"))

    p2 = plot(betas, drift_extinct_probs, ribbon=se_drift_extinct_list, fillalpha=0.3, label="Empirical Drift Extinction", legend=:bottomleft, yscale=:log10, xscale=:log10)
    plot!(p2, betas, theoretical_drift_extinct_probs, line=:dash, color=:blue, label="Theoretical Drift Extinction")
    plot!(p2, betas, immunity_extinct_probs, ribbon=se_immunity_extinct_list, fillalpha=0.3, label="Immunity Extinction")
    plot!(p2, betas, survival_probs, ribbon=se_survival_list, fillalpha=0.3, label="Survival")
    
    # Format the plot
    plot!(p2, xlabel="Beta", ylabel="Probability", title="Outcomes as a function of Beta", size=(800, 600))

    # Save the plot
    savefig(p2, joinpath(OUTPUT_PATH, "outcomes_vs_beta_loglog.png"))
end

main()
