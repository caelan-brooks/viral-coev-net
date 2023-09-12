module CoevolutionNetworkBase
export Population, Network, Simulation, run_simulation!, calculate_total_infected

using Random
using Distributions

"""
    Population

Description:
    A structure representing a population in a simulation, which contains various properties such as 
    viral and immune density distributions, parameters for the population dynamics, and other properties.

Fields:
    - L (Float64): The physical size of the population domain.
    - dx (Float64): The size of a single step in the spatial discretization.
    - r (Float64): A parameter representing ...
    - M (Int): A parameter representing ...
    - beta (Float64): The transmission rate of the virus.
    - alpha (Float64): The recovery rate from the virus.
    - gamma (Float64): The rate at which ...
    - D (Float64): The diffusion coefficient, representing the rate of spatial spread of the virus.
    - Nh (Int): The total number of individuals in the population.
    - viral_density (Vector{Float64}): The density of the virus at various spatial locations.
    - immune_density (Vector{Float64}): The density of immune responses at various spatial locations.
    - stochastic (Bool): A flag indicating whether the simulation should include stochastic effects.
    - time_stamp (Float64): The current time in the simulation.
    - xs (Vector{Float64}): A vector representing the spatial discretization points.
    - num_antigen_points (Int): The number of points in the antigen grid.

Usage:
    pop = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density; stochastic=true, time_stamp=0.0)
"""
struct Population
    L::Float64
    dx::Float64
    r::Float64
    M::Int
    beta::Float64
    alpha::Float64
    gamma::Float64
    D::Float64
    Nh::Int
    viral_density::Vector{Float64}
    immune_density::Vector{Float64}
    stochastic::Bool
    time_stamp::Float64
    xs::Vector{Float64}
    num_antigen_points::Int
end

"""
    Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density; stochastic=true, time_stamp=0.0)

Description:
    Constructor function for creating an instance of the Population structure. 

Parameters:
    - L (Float64): The physical size of the population domain.
    - dx (Float64): The size of a single step in the spatial discretization.
    ... (similarly describe other parameters)

Returns:
    An instance of the Population structure with the given parameters and initialized values for xs and num_antigen_points.

Examples:
    pop = Population(1.0, 0.1, 0.5, 10, 0.3, 0.2, 0.1, 0.05, 100, [0.0 for i in 1:10], [0.0 for i in 1:10])

"""
function Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density; stochastic=true, time_stamp=0.0)
    xs = collect(-L/2:dx:L/2-dx)  # Creating a vector of spatial discretization points
    num_antigen_points = length(xs)  # Calculating the number of points in the antigen grid
    return Population(L, dx, r, M, beta, alpha, gamma, D, Nh, copy(viral_density), copy(immune_density), stochastic, time_stamp, xs, num_antigen_points)
end



"""
    Network

Description:
    A structure representing a network of populations. This is used to simulate the interactions
    between different populations, including migration behaviors which are governed by the migration matrix.

Fields:
    - populations (Vector{Population}): A vector containing the population objects representing different groups in the network.
    - migration_matrix (Matrix{Float64}): A matrix representing the migration rates between different populations in the network.

Constructor:
    Network(populations::Vector{Population}, migration_matrix::Matrix{Float64})

Usage:
    net = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])

Raises:
    - DomainError: If the dimensions of the migration matrix do not match the number of populations.

Examples:
    pop1 = Population(1.0, 0.1, 0.5, 10, 0.3, 0.2, 0.1, 0.05, 100, [0.0 for i in 1:10], [0.0 for i in 1:10])
    pop2 = Population(1.0, 0.1, 0.5, 10, 0.3, 0.2, 0.1, 0.05, 100, [0.0 for i in 1:10], [0.0 for i in 1:10])
    net = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])

"""
struct Network
    populations::Vector{Population}
    migration_matrix::Matrix{Float64}

    """
        Network(populations::Vector{Population}, migration_matrix::Matrix{Float64})

    Description:
        Constructor function for creating an instance of the Network structure. It validates that the 
        dimensions of the migration matrix match the number of populations before creating a new instance.

    Parameters:
        - populations (Vector{Population}): A vector of Population instances representing different groups in the network.
        - migration_matrix (Matrix{Float64}): A matrix representing the migration rates between different populations.

    Returns:
        A new instance of the Network structure with the given populations and migration matrix.

    Raises:
        - DomainError: If the dimensions of the migration matrix do not match the number of populations.
    """
    function Network(populations::Vector{Population}, migration_matrix::Matrix{Float64})
        if size(migration_matrix, 1) != length(populations) || size(migration_matrix, 2) != length(populations)
            throw(DomainError("Migration matrix dimensions do not match the number of populations"))
        end
        new(populations, migration_matrix)
    end
end


"""
    Simulation

Description:
    A mutable structure representing a simulation of network evolution over time. 
    It stores the initial state of the network and the parameters for the simulation time-step and duration. 
    The results of the simulation are stored in the trajectory field.

Fields:
    - initial_network (Network): The initial state of the network before the simulation starts.
    - dt (Float64): The time step used for the simulation.
    - duration (Float64): The total duration for which the simulation will run.
    - trajectory (Vector{Network}): A vector to store the evolution of the network at each time step.
    - duration_times (Vector{Float64}): A vector representing the time points at which the network state is recorded.
    - simulation_complete (Bool): A flag indicating whether the simulation has been completed.

Usage:
    sim = Simulation(initial_network, dt, duration, trajectory, duration_times, false)

Examples:
    init_network = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])
    sim = Simulation(init_network, 0.1, 10.0, [init_network], collect(0:0.1:10), false)

"""
mutable struct Simulation
    initial_network::Network  # The initial state of the network
    dt::Float64  # The time step for the simulation
    duration::Float64  # The total duration of the simulation
    trajectory::Vector{Network}  # A vector to store the network state at each time step
    duration_times::Vector{Float64}  # The time points at which the network state is recorded
    simulation_complete::Bool  # A flag to indicate if the simulation is complete
end


"""
    Simulation(initial_network::Network, dt::Float64, duration::Float64)

Description:
    A constructor for the Simulation struct. This constructor initializes a Simulation instance
    with the given initial network, time step, and duration. The initial state of the network is 
    copied into the trajectory, and the simulation_complete flag is set to false, indicating that
    the simulation has not yet been run.

Parameters:
    - initial_network (Network): The initial state of the network before the simulation starts.
    - dt (Float64): The time step used for the simulation.
    - duration (Float64): The total duration for which the simulation will run.

Returns:
    - A new instance of the Simulation struct initialized with the provided parameters and with 
      simulation_complete flag set to false.

Examples:
    init_network = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])
    sim = Simulation(init_network, 0.1, 10.0)

"""
function Simulation(initial_network::Network, dt::Float64, duration::Float64)
    # Create a range representing the time points at which the network state will be recorded
    duration_times = 0:dt:duration

    # Initialize the trajectory with the initial network state
    trajectory = [copy(initial_network)]

    # Create and return a new Simulation instance with the initial network, 
    # time step, duration, trajectory, and duration times
    # The simulation_complete flag is initially set to false
    Simulation(initial_network, dt, duration, trajectory, duration_times, false)
end


"""
    run_simulation!(sim::Simulation)

Description:
    This function runs the simulation on a given `Simulation` instance. If the simulation has already
    been completed (as indicated by the simulation_complete flag), it raises an error to prevent re-running 
    the simulation. Otherwise, it iteratively evolves the network using the `single_step_evolve_network` 
    function at each time step, recording the new network state in the trajectory. Once the simulation 
    completes, it sets the simulation_complete flag to true.

Parameters:
    - sim (Simulation): An instance of the Simulation struct to be run.

Returns:
    - Modifies the input `sim` in-place, updating the trajectory and simulation_complete flag.

Examples:
    # Initializing a Network and Simulation instance first
    init_network = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])
    sim = Simulation(init_network, 0.1, 10.0)
    # Running the simulation
    run_simulation!(sim)

Raises:
    - Error if the simulation has already been completed.

"""
function run_simulation!(sim::Simulation)
    # Check if the simulation has already been completed, to prevent re-running
    if sim.simulation_complete
        error("Simulation has already been run and cannot be run again.")
    end
    
    # Iteratively evolve the network at each time step in the duration_times (skipping the initial time)
    for time in sim.duration_times[2:end]
        # Get the new network state by evolving the current state by one time step
        new_network = single_step_evolve_network(copy(sim.initial_network), sim.dt)
        
        # Add the new network state to the trajectory
        push!(sim.trajectory, new_network)
        
        # Update the current network state for the next iteration
        sim.initial_network = new_network
    end

    # Mark the simulation as completed
    sim.simulation_complete = true
end


"""
    extract_time_stamps(sim::Simulation)

Description:
    This function extracts the time stamps from each network state in the simulation's trajectory.
    It specifically retrieves the time stamps from the first population of each network in the trajectory.

Parameters:
    - sim (Simulation): An instance of the Simulation struct from which the time stamps will be extracted.

Returns:
    - A vector of time stamps retrieved from the first population of each network state in the trajectory.

Examples:
    # Initializing a Network and Simulation instance first
    init_network = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])
    sim = Simulation(init_network, 0.1, 10.0)
    run_simulation!(sim)
    # Extracting the time stamps
    time_stamps = extract_time_stamps(sim)

"""
function extract_time_stamps(sim::Simulation)
    # Iterating over each network in the simulation's trajectory 
    # and extracting the time stamp from the first population of each network
    return [network.populations[1].time_stamp for network in sim.trajectory]
end

"""
    calculate_total_infected(simulation::Simulation)

Description:
    This function calculates the total number of infected individuals at each time step in the simulation's 
    trajectory. It iterates over each network state in the trajectory and sums the product of the 
    viral_density and dx (spacing) for each population in the network.

Parameters:
    - simulation (Simulation): An instance of the Simulation struct, which contains the trajectory of network 
                               states to analyze.

Returns:
    - A vector containing the total number of infected individuals at each time step in the simulation's 
      trajectory.

Examples:
    # Initializing a Network and Simulation instance first
    init_network = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])
    sim = Simulation(init_network, 0.1, 10.0)
    run_simulation!(sim)
    # Calculating the total number of infected individuals
    total_infected = calculate_total_infected(sim)

"""
function calculate_total_infected(simulation::Simulation)
    # Initialize a vector to hold the total number of infected individuals at each time step
    total_infected = zeros(Float64, length(simulation.trajectory))

    # Iterate over each network state in the simulation's trajectory
    for (i, network) in enumerate(simulation.trajectory)
        # Iterate over each population in the current network state
        for population in network.populations
            # Increment the total infected count for the current time step 
            # by adding the integral of the viral density (viral_density multiplied by dx)
            total_infected[i] += sum(population.viral_density .* population.dx)
        end
    end

    # Return the vector containing the total number of infected individuals at each time step
    return total_infected
end


# function single_step_evolve_network(network::Network, dt)
#     new_populations = [single_step_evolve(copy(pop), dt) for pop in network.populations]

#     # Create a list to store the new viral densities after considering migration effects
#     new_viral_densities = [copy(pop.viral_density) for pop in new_populations]

#     for i in 1:length(new_populations)
#         migration_effect = calculate_migration_effect(new_populations, network.migration_matrix, i)
#         new_viral_densities[i] .+= dt .* migration_effect
#     end

#     # Assign the new viral densities back to the populations
#     for i in 1:length(new_populations)
#         new_populations[i].viral_density = new_viral_densities[i]
#     end

#     # Return a new Network instance with the updated populations and the same migration matrix
#     return Network(new_populations, copy(network.migration_matrix))
# end

"""
    single_step_evolve_network(network::Network, dt)

Description:
    This function evolves the network by a single time step. It first evolves each population 
    in the network independently using the `single_step_evolve` function. Then, it considers 
    the effects of migration between populations on the viral densities. A new Network instance 
    is created with the updated populations and the same migration matrix.

Parameters:
    - network (Network): An instance of the Network struct, which contains the populations and 
                         migration matrix to consider.
    - dt (Float64): The time step over which to evolve the network.

Returns:
    - A new Network instance with updated populations based on both the individual evolution 
      of each population and the migration effects between populations.

Examples:
    # Initializing a Network instance first
    init_network = Network([pop1, pop2], [[0.1, 0.2], [0.3, 0.4]])
    # Evolving the network by a single time step
    new_network = single_step_evolve_network(init_network, 0.1)

"""
function single_step_evolve_network(network::Network, dt)
    # Evolving each population independently using the single_step_evolve function
    new_populations = [single_step_evolve(copy(pop), dt) for pop in network.populations]

    # Initializing a list to store the new viral densities considering migration effects
    new_viral_densities = [copy(pop.viral_density) for pop in new_populations]

    # Iterating over each population in the new populations list
    for i in 1:length(new_populations)
        # Calculating the migration effect on the current population
        migration_effect = calculate_migration_effect(new_populations, network.migration_matrix, i)
        
        # Updating the viral density of the current population based on the migration effect
        new_viral_densities[i] .+= dt .* migration_effect

        # Creating a new Population instance with the updated viral density and the same properties 
        # as the current population in new_populations list
        new_populations[i] = Population(
            new_populations[i].L,
            new_populations[i].dx,
            new_populations[i].r,
            new_populations[i].M,
            new_populations[i].beta,
            new_populations[i].alpha,
            new_populations[i].gamma,
            new_populations[i].D,
            new_populations[i].Nh,
            new_viral_densities[i],  # Updated viral density
            copy(new_populations[i].immune_density),  # Keeping the same immune density
            new_populations[i].stochastic,
            new_populations[i].time_stamp,
            new_populations[i].xs,  # Keeping the same antigenic points
            new_populations[i].num_antigen_points
        )
    end

    # Returning a new Network instance with the updated populations and the same migration matrix
    return Network(new_populations, copy(network.migration_matrix))
end


"""
    single_step_evolve(population::Population, dt)

Description:
    This function performs a single step evolution on a population within a network. 
    It computes the change in viral and immune densities over a time step `dt` by considering 
    various biological and immunological factors such as mutation effect, cross-reactive 
    convolutions, susceptibility, fitness, and growth rate. It also considers stochastic 
    variations in viral densities if the `stochastic` attribute of the population is true. 

Parameters:
    - population (Population): An instance of the Population struct, which contains the details 
                               of the population to evolve.
    - dt (Float64): The time step over which to evolve the population.

Returns:
    - A new Population instance with updated viral and immune densities and an incremented 
      time stamp.

Example:
    # Initializing a Population instance first
    init_pop = Population( ... )  # Add appropriate arguments
    # Evolving the population by a single time step
    new_pop = single_step_evolve(init_pop, 0.1)
"""
function single_step_evolve(population::Population, dt)
    viral_density = copy(population.viral_density)
    immune_density = copy(population.immune_density)

    # Computing the change in viral density due to mutation (diffusion)
    dndt_mutation = compute_mutation_effect(population.D, population.dx, viral_density)

    # Computing the cross-reactive field using the convolution method 
    cross_reactive = cross_reactive_convolution(population.num_antigen_points, immune_density, population.dx, population.r)

    # Computing the susceptibility at each antigenic point based on the cross-reactive field
    susceptibility = compute_susceptibility(population.num_antigen_points, cross_reactive, population.M)

    # Computing the fitness at each antigenic point based on the susceptibility
    fitness = compute_fitness(population.beta, susceptibility, population.alpha, population.gamma)

    # Computing the growth rate of the viral population based on the fitness
    dndt_growth = compute_growth_rate(fitness, viral_density)

    # Computing the total viral population size
    total_viral_pop = compute_total_viral_pop(population.dx, viral_density)

    # Computing the change in immune density based on the current viral and immune densities
    dhdt = compute_immune_density_change(viral_density, total_viral_pop, immune_density, population.M, population.Nh)

    # Updating the viral density using the Euler method
    viral_density .+= dt .* (dndt_mutation .+ dndt_growth)

    # Updating the immune density using the Euler method
    immune_density .+= dt .* dhdt

    # Applying stochastic variations to the viral densities if the `stochastic` parameter is True
    if population.stochastic
        viral_density[viral_density .< 0] .= 0  # Ensuring viral densities remain non-negative
        viral_density = apply_stochasticity(population.dx, viral_density)
    end

    # Creating a new Population instance with updated attributes
    new_population = Population(
        population.L, 
        population.dx, 
        population.r, 
        population.M, 
        population.beta, 
        population.alpha, 
        population.gamma, 
        population.D, 
        population.Nh, 
        viral_density, 
        immune_density, 
        stochastic = population.stochastic, 
        time_stamp = population.time_stamp + dt  # Incrementing the time stamp
    )

    return new_population
end


# Computes the mutation effect based on the given viral density and other parameters.
function compute_mutation_effect(D, dx, viral_density)
    return D / dx^2 * (circshift(viral_density, 1) + circshift(viral_density, -1) - 2 * viral_density)
end

# Computes the susceptibility at each antigenic point based on the cross-reactive convolution and other parameters.
function compute_susceptibility(num_antigen_points, cross_reactive, M)
    return (ones(num_antigen_points) - cross_reactive) .^ M
end

# Computes the fitness at each antigenic point based on susceptibility and other parameters.
function compute_fitness(beta, susceptibility, alpha, gamma)
    return beta .* susceptibility - alpha * ones(length(susceptibility)) - gamma * ones(length(susceptibility))
end

# Computes the growth rate of the viral population based on fitness and viral density.
function compute_growth_rate(fitness, viral_density)
    return fitness .* viral_density
end

# Computes the total size of the viral population by summing up the viral densities.
function compute_total_viral_pop(dx, viral_density)
    return sum(viral_density) * dx
end

# Computes the change in immune density based on viral density, total viral population, and other parameters.
function compute_immune_density_change(viral_density, total_viral_pop, immune_density, M, Nh)
    return 1/(M * Nh) * (viral_density - total_viral_pop .* immune_density)
end

# Applies stochastic variations to the viral densities using a Poisson distribution based on given viral density and other parameters.
function apply_stochasticity(dx, viral_density)
    new_viral_density = rand.(Poisson.(dx .* viral_density)) ./ dx
    return new_viral_density
end


function cross_reactive_convolution(num_antigen_points, immune_density, dx, r)
    """
    Returns the cross-reactive field c(x,t).

    Parameters:
    num_antigen_points (int): Number of antigenic points.
    immune_density (Array): The immune density at each antigenic point.
    dx (float): Discretization of antigenic space.
    r (float): Cross-reactivity parameter.

    Returns:
    Array: The cross-reactive field at each antigenic point.
    """
    
    cross_reactive = zeros(num_antigen_points)  # Initialize an array with zeros to hold the cross-reactive values
    for i in 1:num_antigen_points  # Loop through each antigenic point (1-indexed in Julia)
        for j in 1:num_antigen_points  # For each antigenic point, loop through all other antigenic points
            # Calculate the minimum distance between the current pair of antigenic points, considering the periodic boundary conditions
            diff = min(abs(i - j), num_antigen_points - abs(i - j)) * dx  
            # Increment the cross-reactive value for the i-th point based on the contribution from the j-th point
            cross_reactive[i] += immune_density[j] * exp(-diff/r) * dx  
        end
    end

    return cross_reactive  # Return the computed cross-reactive field
end

# Function to validate the dimensions of the migration matrix
function validate_migration_matrix(matrix, populations)
    if size(matrix, 1) != length(populations) || size(matrix, 2) != length(populations)
        error("Migration matrix dimensions do not match the number of populations.")
    end
end

# Function to calculate the effect of migration on the viral density of a specific population
function calculate_migration_effect(populations, migration_matrix, idx)
    migration_effect = zeros(length(populations[idx].viral_density))

    for j in 1:length(populations)
        if idx != j
            migration_rate_in = migration_matrix[idx, j]
            migration_rate_out = migration_matrix[j, idx]
            migration_effect .+= migration_rate_in .* populations[j].viral_density - migration_rate_out .* populations[idx].viral_density
        end
    end

    return migration_effect
end


"""
    copy(population::Population)

Creates a deep copy of the Population instance.

# Parameters
- `population::Population`: The Population instance to copy.

# Returns
- A new Population instance with the same properties as the original instance.
"""
function Base.copy(population::Population)
    return Population(
        copy(population.L),
        copy(population.dx),
        copy(population.r),
        copy(population.M),
        copy(population.beta),
        copy(population.alpha),
        copy(population.gamma),
        copy(population.D),
        copy(population.Nh),
        copy(population.viral_density), # Ensuring the arrays are deeply copied
        copy(population.immune_density), # Ensuring the arrays are deeply copied
        stochastic = copy(population.stochastic),
        time_stamp = copy(population.time_stamp)
    )
end

function Base.copy(network::Network)
    new_populations = [copy(pop) for pop in network.populations]
    return Network(new_populations, copy(network.migration_matrix))
end

end