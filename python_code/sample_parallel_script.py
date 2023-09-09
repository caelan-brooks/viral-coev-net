import numpy as np
from multiprocessing import Pool
from coevolution_network_base import Population, Network, Simulation
import pickle

# Define function to run a single simulation
def run_single_simulation(simulation_number):
    # The parameters and initializations here are the same as in your original script
    L = 40.0
    dx = 0.3
    x = np.arange(-L/2, L/2, dx)
    r = 3
    M = 15
    beta = 2
    alpha = 1
    gamma = 0
    D = 0.01
    Nh = 10**6
    dt = 0.05
    duration = 80

    viral_density = np.where(np.abs(x) <= 0.5 , 100.0, 0)
    viral_density2 = np.zeros_like(viral_density)
    immune_density = np.zeros_like(viral_density)

    population1 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density)
    population2 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density2, immune_density)

    migration_matrix = np.array([[0, 0], 
                                 [0, 0]])

    network = Network([population1, population2], migration_matrix)
    simulation = Simulation(network, dt, duration)
    simulation.run_simulation()

    # Save the entire simulation object using pickle
    with open(f'simulation_results_{simulation_number}.pkl', 'wb') as file:
        pickle.dump(simulation, file)
    
    

if __name__ == '__main__':
    # Creating a pool of workers to run simulations in parallel
    with Pool() as pool:
        pool.map(run_single_simulation, range(30))