import numpy as np
from multiprocessing import Pool
from coevolution_network_base import Population, Network, Simulation, calculate_total_infected
import pickle
import os

if not os.path.exists("simresults_largedt"):
    os.makedirs("simresults_largedt")


# Define function to run a single simulation
def run_single_simulation(args):
    print(args)
    migration_rate, simulation_number = args

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
    
    migration_matrix = np.array([[0, migration_rate], 
                                 [migration_rate, 0]])

    network = Network([population1, population2], migration_matrix)
    simulation = Simulation(network, dt, duration)
    simulation.run_simulation()

    total_infected = calculate_total_infected(simulation)

    result_dict = {'times': simulation.times, 'total_infected_number': total_infected}
    
    with open(f'simresults_largedt/simulation_results_migration_{migration_rate}_replicate_{simulation_number}.pkl', 'wb') as file:
        pickle.dump(result_dict, file)

if __name__ == '__main__':
    # migration_rates = np.logspace(-6,-0.5,8)  # Example migration rates to sweep over
    migration_rates = np.logspace(-6,0.5,9)  # Example migration rates to sweep over
    start_rep = 0
    num_replicates = 1000

    # Creating a list of tuples with migration rates and simulation numbers
    simulation_args = [(rate, num) for rate in migration_rates for num in range(start_rep, start_rep + num_replicates)]

    with Pool() as pool:
        pool.map(run_single_simulation, simulation_args)
