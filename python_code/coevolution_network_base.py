import numpy as np
from numba import jit
from matplotlib import rc
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

class Network:
    def __init__(self, populations, migration_matrix):
        """
        Initialize a network of populations.

        Parameters:
        populations (list): A list of Population objects.
        migration_matrix (np.array): An nxn matrix representing migration rates 
                                     between populations, where n is the number of 
                                     Population objects.
        """
        self.populations = populations
        self.migration_matrix = migration_matrix
        validate_migration_matrix(migration_matrix, populations)
    
    def single_step_evolve_network(self, dt):
        """
        Evolves each population by a single time step and incorporates migration effects.

        Parameters:
        dt (float): The time step size.
        """
        
        # First, evolve each population by a single time step using their single_step_evolve method
        for pop in self.populations:
            pop.single_step_evolve(dt)
        
        # Then, calculate the change in viral densities due to migration
        new_viral_densities = [pop.viral_density for pop in self.populations]  # Create a copy to store new densities
        for i, pop in enumerate(self.populations):
            migration_effect = calculate_migration_effect(self.populations, self.migration_matrix, i)
            new_viral_densities[i] += dt * migration_effect
        
        # Assign the new viral densities back to the populations
        for i, pop in enumerate(self.populations):
            pop.viral_density = new_viral_densities[i]
    



class Population:
    '''
    This class holds all the data, viral densities and immune densities for a given deme. We will use this class to construct
    a Network object as a collection of Population objects
    '''
    def __init__(self, L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density, stochastic=True, time_stamp=0.0):
        self.L = L  # Length of antigenic space
        self.dx = dx  # Discretization of antigenic space
        self.r = r  # Cross-reactivity
        self.M = M  # Number of immune pressures per host
        self.beta = beta  # Transmission rate
        self.alpha = alpha  # Death rate
        self.gamma = gamma  # Recover rate
        self.D = D  # Mutation rate (diffusion rate)
        self.Nh = Nh  # Number of hosts in this population
        self.viral_density = viral_density.copy() # viral density in this population at t = 0
        self.immune_density = immune_density.copy() # immune density in this population at t=0
        self.stochastic = stochastic # is the simulation stochastic? 
        self.time_stamp = time_stamp

        self.xs = np.arange(-L/2, L/2, dx)  # Vector of antigenic points
        self.num_antigen_points = np.size(self.xs)  # Number of antigenic points
        

    def single_step_evolve(self, dt):
        self.viral_density, self.immune_density = single_step_evolve(
            dt, 
            D=self.D, 
            dx=self.dx, 
            viral_density=self.viral_density, 
            beta=self.beta, 
            alpha=self.alpha, 
            gamma=self.gamma, 
            M=self.M, 
            num_antigen_points=self.num_antigen_points, 
            immune_density=self.immune_density, 
            Nh=self.Nh, 
            r=self.r, 
            cross_reactive_convolution_func=cross_reactive_convolution,
            stochastic = self.stochastic
        )

        self.time_stamp += dt

    def plot(self, ax1=None, ax2=None, color1='tab:blue', color2='tab:red'):
        if ax1 is None and ax2 is None:
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()

        ax1.set_xlabel('Antigenic coordinate $x$')
        ax1.set_ylabel('Viral Density', color=color1)
        ax1.set_xlim(min(self.xs), max(self.xs))
        ax1.set_ylim(1, self.Nh)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1.set_yscale('log')
        ax1.plot(self.xs, self.viral_density, color=color1, label='Viral Density')
        ax1.legend(loc='upper left')

        ax2.set_ylabel('Immune Density', color=color2)
        ax2.set_ylim(0, 1)
        ax2.tick_params(axis='y', labelcolor=color2)
        ax2.plot(self.xs, self.immune_density, color=color2, label='Immune Density')
        ax2.legend(loc='upper right')

        plt.show()

@jit(nopython=True)
def cross_reactive_convolution(num_antigen_points, immune_density, dx, r):
    """
    Returns the cross-reactive field c(x,t).

    Parameters:
    num_antigen_points (int): Number of antigenic points.
    immune_density (np.array): The immune density at each antigenic point.
    dx (float): Discretization of antigenic space.
    r (float): Cross-reactivity parameter.

    Returns:
    np.array: The cross-reactive field at each antigenic point.
    """
    
    cross_reactive = np.zeros(num_antigen_points)  # Initialize an array with zeros to hold the cross-reactive values
    for i in range(num_antigen_points):  # Loop through each antigenic point
        for j in range(num_antigen_points):  # For each antigenic point, loop through all other antigenic points
            # Calculate the minimum distance between the current pair of antigenic points, considering the periodic boundary conditions
            diff = min(np.abs(i - j), num_antigen_points - np.abs(i - j)) * dx  
            # Increment the cross-reactive value for the i-th point based on the contribution from the j-th point
            cross_reactive[i] += immune_density[j] * np.exp(-diff/r) * dx  

    return cross_reactive  # Return the computed cross-reactive field

@jit(nopython=True)
def single_step_evolve(dt, D, dx, viral_density, beta, alpha, gamma, M, num_antigen_points, immune_density, Nh, r, cross_reactive_convolution_func, stochastic=True):
    """
    Evolves the viral_density and immune_density fields by a single time-step dt using the Euler method.

    Parameters:
    dt (float): The time step size.
    D (float): Mutation rate (diffusion rate).
    dx (float): Discretization of antigenic space.
    viral_density (np.array): The viral density at each antigenic point at time t.
    beta (float): Transmission rate.
    alpha (float): Death rate.
    gamma (float): Recover rate.
    M (int): Number of immune pressures per host.
    num_antigen_points (int): Number of antigenic points.
    immune_density (np.array): The immune density at each antigenic point at time t.
    Nh (int): Number of hosts in the population.
    r (float): Cross-reactivity parameter.
    cross_reactive_convolution_func (function): The function to compute the cross-reactive convolution.

    Returns:
    tuple: The updated viral and immune densities.
    """
    
    # Compute the change in viral density due to mutation (diffusion)
    dndt_mutation = compute_mutation_effect(D, dx, viral_density)
    
    # Compute the cross-reactive field using the provided function
    cross_reactive = cross_reactive_convolution_func(num_antigen_points, immune_density, dx, r)
    
    # Compute the susceptibility at each antigenic point based on the cross-reactive field
    susceptibility = compute_susceptibility(num_antigen_points, cross_reactive, M)
    
    # Compute the fitness at each antigenic point based on the susceptibility
    fitness = compute_fitness(beta, susceptibility, alpha, gamma)
    
    # Compute the growth rate of the viral population based on the fitness
    dndt_growth = compute_growth_rate(fitness, viral_density)
    
    # Compute the total viral population size
    total_viral_pop = compute_total_viral_pop(dx, viral_density)
    
    # Compute the change in immune density based on the current viral and immune densities
    dhdt = compute_immune_density_change(viral_density, total_viral_pop, immune_density, M, Nh)
    
    # Update the viral density using the Euler method
    viral_density += dt * (dndt_mutation + dndt_growth)
    
    # Update the immune density using the Euler method
    immune_density += dt * dhdt
    
    # Apply stochastic variations to the viral densities if stochastic parameter is True
    if stochastic:
        viral_density = apply_stochasticity(dx, num_antigen_points, viral_density)
    
    return viral_density, immune_density  # Return the updated viral and immune densities

@jit(nopython=True)
def compute_mutation_effect(D, dx, viral_density):
    return D / dx**2 * (np.roll(viral_density, 1) + np.roll(viral_density, -1) - 2 * viral_density)

@jit(nopython=True)
def compute_susceptibility(num_antigen_points, cross_reactive, M):
    return np.power(np.ones(num_antigen_points) - cross_reactive, M)

@jit(nopython=True)
def compute_fitness(beta, susceptibility, alpha, gamma):
    return beta * susceptibility - alpha - gamma

@jit(nopython=True)
def compute_growth_rate(fitness, viral_density):
    return fitness * viral_density

@jit(nopython=True)
def compute_total_viral_pop(dx, viral_density):
    return np.sum(viral_density) * dx

@jit(nopython=True)
def compute_immune_density_change(viral_density, total_viral_pop, immune_density, M, Nh):
    return 1/(M * Nh) * (viral_density - total_viral_pop * immune_density)

@jit(nopython=True)
def apply_stochasticity(dx, num_antigen_points, viral_density):
    for i in range(num_antigen_points):
        viral_density[i] = np.random.poisson(dx * viral_density[i]) / dx
    return viral_density

# Utility Functions for migration calculation
def validate_migration_matrix(matrix, populations):
    if matrix.shape[0] != len(populations) or matrix.shape[1] != len(populations):
        raise ValueError("Migration matrix dimensions do not match the number of populations.")

def calculate_migration_effect(populations, migration_matrix, idx):
    migration_effect = np.zeros_like(populations[idx].viral_density)
    for j, other_pop in enumerate(populations):
        if idx != j:
            migration_rate = migration_matrix[idx, j]
            migration_effect += migration_rate * other_pop.viral_density
    return migration_effect


