module MultistrainSIR
export ModelParameters, initialize_state, simulate_trajectories

using Random
using Distributions

# Define the parameters for your model
struct ModelParameters
    β::Matrix{Float64}  # Transmission rates for each strain in each deme
    γ::Matrix{Float64}  # Recovery rates for each strain in each deme
    k::Matrix{Float64}  # Migration matrix between demes
    μ::Matrix{Float64}  # Mutation matrix between strains
    σ::Matrix{Float64}  # Scaling factor for demographic noise
    N_h::Vector{Int}    # Total host population in each deme
    dt::Float64         # Time step for the simulation
    total_time::Float64 # Total time for the simulation
end


# Define a function to initialize the state variables
function initialize_state(num_demes::Int, num_strains::Int)
    S = ones(num_demes, num_strains)  # Susceptibility initialized to 1 for all
    N = zeros(num_demes, num_strains) # Number of infected individuals
    return S, N
end

using Random  # Make sure to include this at the beginning of your module for random number generation

function euler_step(S::Matrix{Float64}, N::Matrix{Float64}, params::ModelParameters)
    # Create copies of S and N to avoid altering the original ones
    S_new = copy(S)
    N_new = copy(N)
    
    # Unpack parameters
    β, γ, k, μ, σ, N_h = params.β, params.γ, params.k, params.μ, params.σ, params.N_h
    dt = params.dt

    num_demes, num_strains = size(β)
    dS = zeros(num_demes, num_strains)
    dN = zeros(num_demes, num_strains)

    # Compute the changes using the SIR model dynamics
    for i in 1:num_demes
        for j in 1:num_strains
            # Susceptible dynamics
            dS[i, j] = -β[i, j] / N_h[i] * N_new[i, j] * S_new[i, j] * dt

            # Infected dynamics
            infection = (β[i, j] * S_new[i, j] * N_new[i, j] - γ[i, j] * N_new[i, j]) * dt

            # Migration dynamics
            migration_in = sum(k[i, :] .* N_new[:, j]) * dt
            migration_out = sum(k[:, i] .* N_new[i, j]) * dt
            migration = migration_in - migration_out

            # Mutation dynamics
            mutation_in = sum(μ[j, :] .* N_new[i, :]) * dt
            mutation_out = sum(μ[:, j] .* N_new[i, j]) * dt
            mutation = mutation_in - mutation_out

            # Net change in the number of infected individuals with strain j in deme i
            dN[i, j] = infection + migration + mutation
        end
    end

    # Update the state without altering the original S and N
    S_new += dS
    N_new += dN

    # Incorporate demographic noise by redrawing N from a Poisson distribution
    for i in 1:num_demes
        for j in 1:num_strains
            if N_new[i, j] > 0 && σ[i,j] > 0  # Only redraw if there are infected individuals
                λ = max(N_new[i, j] / (σ[i,j]^2 * dt), 0)  # Rate parameter for Poisson distribution
                N_new[i, j] = σ[i,j]^2 * dt * rand(Poisson(λ))
            end
        end
    end

    return S_new, N_new
end

function simulate_trajectories(initial_S::Matrix{Float64}, initial_N::Matrix{Float64}, params::ModelParameters)
    # Determine the number of steps from the total simulation time and dt
    times = 0: params.dt : params.total_time
    num_steps = length(times)

    # Get the number of demes and strains
    num_demes, num_strains = size(initial_S)

    # Initialize the arrays to store the trajectories
    S_trajectory = Array{Float64, 3}(undef, num_demes, num_strains, num_steps)
    N_trajectory = Array{Float64, 3}(undef, num_demes, num_strains, num_steps)

    # Set the initial conditions
    S_current = copy(initial_S)
    N_current = copy(initial_N)

    # Store the initial conditions as the first step in the trajectories
    S_trajectory[:, :, 1] = S_current
    N_trajectory[:, :, 1] = N_current

    # Iterate over each time step
    for step in 2:num_steps
        # Perform the Euler step to get the new S and N
        S_current, N_current = euler_step(S_current, N_current, params)

        # Store the results in the trajectories
        S_trajectory[:, :, step] = S_current
        N_trajectory[:, :, step] = N_current
    end

    return S_trajectory, N_trajectory, times
end

end # module
