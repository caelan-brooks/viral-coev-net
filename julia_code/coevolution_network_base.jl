module CoevolutionNetworkBase
export Population, Network, Simulation, run_simulation!, calculate_total_infected, calculate_total_infected_per_deme, single_step_evolve!, thin_simulation!, plot_spacetime_density, calculate_antigenic_variance_per_deme

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
    sigma::Float64
    Nh::Int
    viral_density::Vector{Float64}
    immune_density::Vector{Float64}
    stochastic::Bool
    time_stamp::Float64
    xs::Vector{Float64}
    num_antigen_points::Int
    temporary_data::Vector{Float64}
    cross_reactive::Vector{Float64}
    susceptibility::Vector{Float64}
    fitness::Vector{Float64}
    noise_method::Symbol
    noise_scaling_factor::Float64
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
function Population(L::Float64, dx::Float64, r::Float64, M::Int64, beta::Float64, alpha::Float64, gamma::Float64, D::Float64, Nh::Int64, viral_density::Vector{Float64}, immune_density::Vector{Float64}; stochastic::Bool=true, time_stamp::Float64=0.0, sigma::Float64=1.0, noise_method::Symbol=:PL_with_dx)
    xs = collect(-L/2:dx:L/2-dx)  # Creating a vector of spatial discretization points
    num_antigen_points = length(xs)  # Calculating the number of points in the antigen grid
    temporary_data = zeros(num_antigen_points)
    cross_reactive = zeros(num_antigen_points)
    susceptibility = zeros(num_antigen_points)
    fitness = zeros(num_antigen_points)

    @assert length(viral_density) == num_antigen_points "Viral density vector size mismatch"
    @assert length(immune_density) == num_antigen_points "Immune density vector size mismatch"

    noise_scaling_factor = (noise_method==:PL) ? 2 / sigma^2 : 2 * dx / sigma^2
    if noise_method == :PL_with_dx
        noise_scaling_factor = 2 * dx / sigma^2
    elseif noise_method == :PL
        noise_scaling_factor = 2 / sigma^2
    else
        error("noise method chosen is not one of the known methods :PL_with_dx or :PL")
    end

    return Population(L, dx, r, M, beta, alpha, gamma, D, sigma, Nh, copy(viral_density), copy(immune_density), stochastic, time_stamp, xs, num_antigen_points, temporary_data,cross_reactive,susceptibility,fitness, noise_method, noise_scaling_factor)
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
    thin_by::Int64 # Number which knows how much time resultion is requested
    cross_reactive_kernel::Matrix{Float64}
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
function Simulation(initial_network::Network, dt::Float64, duration::Float64; thin_by=1)
    # Create a range representing the time points at which the network state will be recorded
    duration_times = collect(0:dt:duration)
    num_time_steps = length(duration_times)

    # Initialize the trajectory with the initial network state
    trajectory = [deepcopy(initial_network) for i in 1:thin_by:num_time_steps]

    # Fill out the cross_reactive_kernel
    dx = initial_network.populations[1].dx
    r = initial_network.populations[1].r
    num_antigen_points = initial_network.populations[1].num_antigen_points
    cross_reactive_kernel = Matrix{Float64}(undef,  num_antigen_points, num_antigen_points)

    for i = 1:num_antigen_points
        for j = 1:num_antigen_points
            diff = min(abs(i - j), num_antigen_points - abs(i - j)) * dx
            cross_reactive_kernel[j,i] = exp(-diff/r) * dx
        end
    end 

    # Create and return a new Simulation instance with the initial network, 
    # time step, duration, trajectory, and duration times
    # The simulation_complete flag is initially set to false
    return Simulation(initial_network, dt, duration, trajectory, duration_times, false, thin_by, cross_reactive_kernel)
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
    
    current_update = 1

    # Iteratively evolve the network at each time step in the duration_times (skipping the initial time)
    for i in 2:length(sim.duration_times)
        single_step_evolve_network!(sim.initial_network, sim.dt, sim.cross_reactive_kernel)
        
        if mod(i-1,sim.thin_by)==0
            current_update += 1
            copy_network_data!(sim.trajectory[current_update],sim.initial_network)
        end
    end

    # Mark the simulation as completed
    sim.simulation_complete = true

    # Make the time vector line up with the saved configurations
    sim.duration_times = sim.duration_times[1:sim.thin_by:end];
end

function copy_network_data!(dest::Network, source::Network)
    # Implement the necessary logic to copy data from the source network to the destination network
    # For example:
    for i in 1:length(dest.populations)
        dest.populations[i].viral_density .= source.populations[i].viral_density
        dest.populations[i].immune_density .= source.populations[i].immune_density
        dest.populations[i].cross_reactive .= source.populations[i].cross_reactive
        dest.populations[i].susceptibility .= source.populations[i].susceptibility
        dest.populations[i].fitness .= source.populations[i].fitness
        # Copy other fields as necessary
    end
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

function calculate_total_infected_per_deme(simulation::Simulation)
    # Get the number of populations (demes)
    num_demes = length(simulation.trajectory[1].populations)
    num_time_points = length(simulation.duration_times)
    
    # Initialize a list of vectors to hold the total number of infected individuals at each time step for each deme
    total_infected_per_deme = zeros(num_demes, num_time_points);

    # Iterate over each network state in the simulation's trajectory
    for (i, network) in enumerate(simulation.trajectory)
        # Iterate over each population in the current network state
        for (j, population) in enumerate(network.populations)
            # Increment the total infected count for the current time step and deme
            # by adding the integral of the viral density (viral_density multiplied by dx)
            total_infected_per_deme[j,i] += sum(population.viral_density .* population.dx)
        end
    end

    # Return the list containing the total number of infected individuals at each time step for each deme
    return total_infected_per_deme
end

function calculate_antigenic_variance_per_deme(simulation::Simulation)
    xs = simulation.trajectory[1].populations[1].xs
    num_time_points = length(simulation.duration_times)
    num_demes = length(simulation.trajectory[1].populations)

    total_infected_per_deme = calculate_total_infected_per_deme(simulation)
    variances_per_deme = zeros(num_demes,num_time_points)

    for i = 1:num_time_points
        for j = 1:num_demes
            population = simulation.trajectory[i].populations[j]
            if total_infected_per_deme[j,i] > 1
                avg_antigenicity = sum(xs .* population.viral_density .* population.dx) ./ total_infected_per_deme[j,i]
                variances_per_deme[j,i] = sum((xs .- avg_antigenicity).^2 .* population.viral_density .* population.dx) ./ total_infected_per_deme[j,i]
            else
                variances_per_deme[j,i] = 0
            end
        
        end
    end

    return variances_per_deme
end
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
function single_step_evolve_network!(network::Network, dt::Float64, cross_reactive_kernel::Matrix{Float64})
    # Evolving each population independently using the single_step_evolve function
    for pop in network.populations
        single_step_evolve!(pop,dt,cross_reactive_kernel)
    end

    calculate_migration_effect!(network,dt)
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
function single_step_evolve!(population::Population, dt::Float64, cross_reactive_kernel::Matrix{Float64})
    # Computing the change in viral density due to mutation (diffusion)
    compute_mutation_effect!(population, dt)

    # Computing the cross-reactive field using the convolution method 
    cross_reactive_convolution!(population, cross_reactive_kernel)

    # Computing the susceptibility at each antigenic point based on the cross-reactive field
    compute_susceptibility!(population)

    # Computing the fitness at each antigenic point based on the susceptibility
    compute_fitness!(population)

    # Updating the immune density using the Euler method
    # total_infected_individuals = sum(population.viral_density) * population.dx
    population.immune_density .+= dt / population.M / population.Nh .* (population.viral_density .- (sum(population.viral_density) * population.dx) .* population.immune_density)
    population.viral_density .+= population.fitness .* population.viral_density .* dt

    # Applying stochastic variations to the viral densities if the `stochastic` parameter is True
    if population.stochastic
        for i in eachindex(population.viral_density)
            if population.viral_density[i] < 0.0
                population.viral_density[i] = 0.0
            end
        end
        
        apply_stochasticity!(population, dt)
    end

end


# Computes the mutation effect based on the given viral density and other parameters.
function compute_mutation_effect!(population::Population, dt::Float64)
    D_over_dx2 = population.D / population.dx^2

    for i in 2:(population.num_antigen_points-1)
        population.temporary_data[i] = D_over_dx2 * dt * (population.viral_density[i-1] + population.viral_density[i+1] - 2 * population.viral_density[i])
    end
    population.temporary_data[1] = D_over_dx2 * dt * (population.viral_density[2] + population.viral_density[end] - 2 * population.viral_density[1])
    population.temporary_data[end] = D_over_dx2 * dt * (population.viral_density[1] + population.viral_density[end-1] - 2 * population.viral_density[end])
    
    population.viral_density .+= population.temporary_data
end


# Computes the susceptibility at each antigenic point based on the cross-reactive convolution and other parameters.
function compute_susceptibility!(population::Population)
    population.susceptibility .= (1.0 .- population.cross_reactive) .^ population.M
end


# Computes the fitness at each antigenic point based on susceptibility and other parameters.
function compute_fitness!(population::Population)
    population.fitness .= population.beta .* population.susceptibility .- population.alpha .- population.gamma
end

# check if this should have dx in it, it seems like the answer is no
function apply_stochasticity!(population::Population, dt::Float64)
    scaling_factor = population.noise_scaling_factor / dt

    for i in eachindex(population.viral_density)
        population.viral_density[i] = rand(Poisson(population.viral_density[i] * scaling_factor)) 
        population.viral_density[i] = (population.viral_density[i] == 0) ? 0 : rand(Gamma(population.viral_density[i])) / scaling_factor 
    end
end



function cross_reactive_convolution!(population::Population, cross_reactive_kernel::Matrix{Float64})
    """
    Modifies the cross-reactive field c(x,t) in place.

    Parameters:
    population (Population): The population object containing the necessary parameters and fields.
    """
    
    population.cross_reactive .= cross_reactive_kernel * population.immune_density
    # population.cross_reactive .= 0.0  # Reset the array to zero without creating a new array
    # for i in 1:population.num_antigen_points  # Loop through each antigenic point (1-indexed in Julia)
    #     for j in 1:population.num_antigen_points  # For each antigenic point, loop through all other antigenic points
    #         # Calculate the minimum distance between the current pair of antigenic points, considering the periodic boundary conditions
    #         diff = min(abs(i - j), population.num_antigen_points - abs(i - j)) * population.dx  
    #         # Increment the cross-reactive value for the i-th point based on the contribution from the j-th point
    #         population.cross_reactive[i] += population.immune_density[j] * exp(-diff/population.r) * population.dx
    #     end
    # end
end



# Function to validate the dimensions of the migration matrix
function validate_migration_matrix(matrix, populations)
    if size(matrix, 1) != length(populations) || size(matrix, 2) != length(populations)
        error("Migration matrix dimensions do not match the number of populations.")
    end
end

function calculate_migration_effect!(network::Network, dt::Float64)
    populations = network.populations
    migration_matrix = network.migration_matrix

    # First, we calculate the migration effects and store them in temporary_data for each population
    for i in 1:length(populations)
        # Resetting the temporary_data before calculating the migration effects
        populations[i].temporary_data .= 0.0

        for j in 1:length(populations)
            if i != j
                # migration_rate_in = migration_matrix[i, j]
                # migration_rate_out = migration_matrix[j, i]
                populations[i].temporary_data .+= migration_matrix[i, j] .* populations[j].viral_density .- migration_matrix[j,i] .* populations[i].viral_density
            end
        end
    end

    # Then, we update the viral_density with the values stored in temporary_data and scale by dt
    for i in 1:length(populations)
        populations[i].viral_density .+= dt .* populations[i].temporary_data
    end
end


function thin_simulation!(sim::Simulation, thin_by::Int64)
    sim.trajectory = [sim.trajectory[1:thin_by:end-1]; sim.trajectory[end]]
    sim.duration_times = [sim.duration_times[1:thin_by:end-1]; sim.duration_times[end]]
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