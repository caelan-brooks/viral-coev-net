using Plots
using Statistics

# Assuming MultistrainSIR is the name of your module and it's in the same directory or installed
include("multistrain_SIR.jl")  # Include the module file if it's in the same directory
using .MultistrainSIR

# Define model parameters
num_demes = 2
num_strains = 2
β = ones(num_demes, num_strains) * 2.5  # Just an example value, adjust accordingly
γ = ones(num_demes, num_strains) * 1.2  # Same here
k = [0.0 1.0; 1.0 0.0] * 1 # Migration rates between the two demes
μ = [0.0 0.0; 1.0 0.0] * 1e-6 # Mutation rates between the two strains
σ = ones(num_demes, num_strains) * 4.5 # Demographic noise intensity
N_h = [10^6, 10^6]  # Total host population in each deme
dt = 0.1  # Time step for Euler integration
total_time = 100.0  # Total simulation time

model_params = ModelParameters(β, γ, k, μ, σ, N_h, dt, total_time)

# Initialize the state
initial_S, initial_N = initialize_state(num_demes, num_strains)
initial_N[1, 1] = 10  # Initial infected in deme 1, strain 1

# Number of simulations
num_simulations = 100

# Run the simulations
@time simulations = [simulate_trajectories(initial_S, initial_N, model_params) for _ in 1:num_simulations]

# Define an array of line styles for the strains
line_styles = [:solid, :dash]

# Initialize an array to store the plots for each deme
deme_plots = []

# Iterate over each deme for plotting
for i in 1:num_demes
    # Initialize a plot for the current deme
    deme_plot = plot(title = "Deme $i", legend = :outertopright, yscale=:log10)

    # Extract and plot the data for each strain within the deme
    for j in 1:num_strains
        # Plot for each simulation run
        for sim in 1:num_simulations
            # Extract the time points, S and N trajectories for this simulation
            S_trajectory, N_trajectory, times = simulations[sim]
            
            # Extract the data for the current strain and deme from the N trajectories
            strain_data = N_trajectory[i, j, :]
            reg = strain_data .> 0

            # Plot the strain data onto the deme plot
            plot!(deme_plot, times[reg], strain_data[reg], label = "Strain $j Run $sim",
                   line = (line_styles[j], 2))  # Use the line style corresponding to the strain
        end
    end

    # Append the deme plot to the array of deme plots
    push!(deme_plots, deme_plot)
end

# Combine all deme plots into a single plot with a layout
final_plot = plot(deme_plots..., layout = (num_demes, 1), size = (600, 400 * num_demes))

# Display the final plot
display(final_plot)
