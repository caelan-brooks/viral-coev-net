using Serialization
using Base.Threads

# Include the MultistrainSIR module
include("multistrain_SIR.jl")
using .MultistrainSIR

# Define the range for migration and mutation rates to sweep over
migration_rates = exp10.(LinRange(-5, 0, 10))
mutation_rates = exp10.(LinRange(-6, -5, 3))
num_simulations = 1 * 10^4    # Number of replicate simulations

# Define the model parameters that remain constant
num_demes = 2
num_strains = 2
β = ones(num_demes, num_strains) * 2.5
γ = ones(num_demes, num_strains) * 1.3
σ = ones(num_demes, num_strains) * 5.0
N_h = [10^6, 10^6]
dt = 0.15 
total_time = 100.0

# Directory for saving the results
results_dir = "simulation_results"
isdir(results_dir) || mkdir(results_dir)

# Function to run simulations and serialize the results
function run_and_serialize_simulations(μ_rate, k_rate, num_simulations, results_dir, k_index, μ_index)
    model_params = ModelParameters(β, γ, [0.0 k_rate; k_rate 0.0], [0.0 μ_rate; μ_rate 0.0], σ, N_h, dt, total_time)
    initial_S, initial_N = initialize_state(num_demes, num_strains)
    initial_N[1, 1] = 10

    # Prepare the filename
    filename = joinpath(results_dir, "trajectories_k$(k_index)_μ$(μ_index).dat")

    # Simulate and collect trajectories
    trajectories = Vector{Any}(undef, num_simulations)
    
    @threads for sim in 1:num_simulations
        _, N_trajectory, _ = simulate_trajectories(initial_S, initial_N, model_params)
        trajectories[sim] = N_trajectory
    end

    # Serialize the trajectories to a file
    open(filename, "w") do io
        serialize(io, trajectories)
    end
end

# Loop over all combinations of migration and mutation rates
for k_index in 1:length(migration_rates)
    for μ_index in 1:length(mutation_rates)
        k_rate = migration_rates[k_index]
        μ_rate = mutation_rates[μ_index]
        
        println("Starting simulations for k_rate=$k_rate, μ_rate=$μ_rate")
        run_and_serialize_simulations(μ_rate, k_rate, num_simulations, results_dir, k_index, μ_index)
        
        println("Finished simulations for k_rate=$k_rate, μ_rate=$μ_rate")
    end
end
