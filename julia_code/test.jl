using Plots
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase
using LaTeXStrings

# Parameters
L = 40.0
dx = 0.3
x = -L/2:dx:L/2-dx
r = 3.0
M = 15
beta = 2.0
alpha = 1.0
gamma = 0.0
D = 0.01
Nh = 10^6

# Initialize viral and immune densities
viral_density = [abs(val) <= 0.5 ? 100.0 : 0.0 for val in x]
viral_density2 = zeros(Float64, length(x))
immune_density = zeros(Float64, length(x))

# Create Population instances
population1 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density)
population2 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density2, immune_density)

# Create Network instance
migration_matrix = [0.0 0.1; 0.1 0.0] # Define an appropriate migration matrix
network = Network([population1, population2], migration_matrix)

# Create Simulation instance
dt = 0.05 # Define an appropriate time step size
duration = 80.0 # Define an appropriate simulation duration
simulation = Simulation(network, dt, duration)

# Define a function to calculate the total number of infected individuals
function calculate_total_infected(simulation)
    total_infected = zeros(Float64, length(simulation.trajectory))
    for (i, network) in enumerate(simulation.trajectory)
        for population in network.populations
            total_infected[i] += sum(population.viral_density .* population.dx)
        end
    end
    return total_infected
end

# Run the simulation
run_simulation!(simulation)

# Calculate total number of infected individuals
total_infected = calculate_total_infected(simulation)

# Plot the results
plotvar = plot(simulation.duration_times, total_infected, 
    label="Total Infected", 
    xlabel="Time", 
    ylabel="Number of Infected Individuals", 
    title="Infection Over Time",
    guidefont=font(16, "Computer Modern"), 
    tickfont=font(14, "Computer Modern"), 
    legendfont=font(14, "Computer Modern"), 
    titlefont=font(16, "Computer Modern"))
display(plotvar)
