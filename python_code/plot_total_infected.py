import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob

def calculate_total_infected(simulation):
    total_infected = np.zeros(len(simulation.trajectory))
    for i, network in enumerate(simulation.trajectory):
        for population in network.populations:
            total_infected[i] += np.sum(population.viral_density * population.dx)
    return total_infected

if __name__ == '__main__':
    # Get a list of all simulation file paths
    simulation_files = glob.glob('simulation_*.pkl')

    # Create a new figure for the plot
    plt.figure()

    count = 0
    # Loop over each file and calculate and plot the total infected number vs time
    for simulation_file in simulation_files:
        count += 1
        print(count)
        # Load the simulation object from the file
        with open(simulation_file, 'rb') as file:
            simulation = pickle.load(file)
        
        # Calculate the total infected number for this simulation
        total_infected = calculate_total_infected(simulation)

        # Plot the total infected number vs time for this simulation
        plt.plot(simulation.times, total_infected, linewidth=2) # Making lines thicker

        # Add a break (NaN value) to prevent connecting lines between different simulations
        plt.plot([np.nan], [np.nan])

    # Add labels and a legend to the plot
    plt.xlabel('Time')
    plt.ylabel('Total Infected Number')
    plt.title('Total Infected Number vs Time')
    
    # Setting y-axis to log scale
    plt.yscale('log')

    # Show the plot
    plt.show()

