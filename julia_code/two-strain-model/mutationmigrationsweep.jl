using Plots
using Statistics
using Distributions
using Base.Threads

# Assuming MultistrainSIR is the name of your module and it's in the same directory or installed
include("multistrain_SIR.jl")  # Include the module file if it's in the same directory
using .MultistrainSIR


# Define the range for migration and mutation rates to sweep over
migration_rates = exp10.(LinRange(-5,0,10))
# mutation_rates = exp10.(LinRange(-8,-6,5))
mutation_rates = [4 * 10^-5]
cutoff = 500  # Define the cutoff for strain two
num_simulations = 10^5    # Number of replicate simulations

# Define the model parameters that remain constant
num_demes = 2
num_strains = 2
β = ones(num_demes, num_strains) * 2.5
γ = ones(num_demes, num_strains) * 1.3
σ = ones(num_demes, num_strains) * 5.0
N_h = [10^6, 10^6]
dt = 0.15
total_time = 100.0

# Initialize the probability matrix
prob_matrix = zeros(length(migration_rates), length(mutation_rates))

# Define a function to run simulations and calculate probability
function calculate_probability(μ_rate, k_rate)
    exceed_count = 0
    model_params = ModelParameters(β, γ, [0.0 k_rate; k_rate 0.0], [0.0 0.0; μ_rate 0.0], σ, N_h, dt, total_time)
    initial_S, initial_N = initialize_state(num_demes, num_strains)
    initial_N[1, 1] = 10

    @threads for _ in 1:num_simulations
        _, N_trajectory, _ = simulate_trajectories(initial_S, initial_N, model_params)
        if any(N_trajectory[1, 2, :] + N_trajectory[1, 2, :] .> cutoff)
            exceed_count += 1
        end
    end
    exceed_count / num_simulations
end

# Assuming `migration_rates` and `mutation_rates` are defined as before
num_migration_rates = length(migration_rates)
num_mutation_rates = length(mutation_rates)
num_combinations = num_migration_rates * num_mutation_rates

# Initialize an array to store the probabilities and their standard errors
prob_matrix = Array{Float64}(undef, num_migration_rates, num_mutation_rates)
stderr_matrix = Array{Float64}(undef, num_migration_rates, num_mutation_rates)

# Multithreaded computation of the probabilities
for k_index in 1:num_migration_rates
    for μ_index in 1:num_mutation_rates
        k_rate = migration_rates[k_index]
        μ_rate = mutation_rates[μ_index]
        prob_matrix[k_index, μ_index] = calculate_probability(μ_rate, k_rate)
        # Calculate the standard error for the probability
        stderr_matrix[k_index, μ_index] = sqrt(prob_matrix[k_index, μ_index] * (1 - prob_matrix[k_index, μ_index]) / num_simulations)
        println("Finished calculation for k_rate=$k_rate, μ_rate=$μ_rate")
    end
end
# Plot the probabilities vs. migration rate for each mutation rate with error ribbons
plot(title="Probability vs. Migration Rate for Different Mutation Rates",
     xlabel="Migration Rate", ylabel="Probability", xscale=:log10,
     legend=:bottomright)

for μ_index in 1:num_mutation_rates
    plot!(migration_rates, prob_matrix[:, μ_index], 
          ribbon=stderr_matrix[:, μ_index],
          label="Mutation Rate $(mutation_rates[μ_index])")
end

# Display the final combined plot
final_plot = current()
display(final_plot)

# Save the combined plot to a file
savefig(final_plot, "probability_vs_migration_rate_with_error.png")
