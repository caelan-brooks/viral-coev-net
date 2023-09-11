import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt

DIRECTORY_PATH = 'simresults_largedt'

def calculate_probability_of_survival(migration_rate, cutoff):
    # Get list of all files with the given migration rate
    files = glob.glob(f'{DIRECTORY_PATH}/simulation_results_migration_{migration_rate}_replicate_*.pkl')

    survival_counts = 0
    total_files = len(files)
    print(total_files)
    for file in files:
        with open(file, 'rb') as f:
            result_dict = pickle.load(f)
        
        # Get the last 5% of the data
        data_slice = result_dict['total_infected_number'][int(-0.05*len(result_dict['total_infected_number'])):]
        
        # Check if the time average over the last 5% of the data exceeds the cutoff
        if np.mean(data_slice) > cutoff:
            survival_counts += 1
    
    # Calculate probability of survival and its standard error
    probability_of_survival = survival_counts / total_files
    standard_error = np.sqrt((probability_of_survival * (1 - probability_of_survival)) / total_files)
    
    return probability_of_survival, standard_error

def plot_total_infected_trajectories(migration_rate):
    files = glob.glob(f'{DIRECTORY_PATH}/simulation_results_migration_{migration_rate}_replicate_*.pkl')
    plt.figure()
    
    for file in files:
        with open(file, 'rb') as f:
            result_dict = pickle.load(f)
        
        plt.plot(result_dict['times'], result_dict['total_infected_number'], label='Replicate')
    
    plt.xlabel('Time')
    plt.ylabel('Total Infected Number')
    plt.title(f'Total Infected Number vs Time (Migration Rate: {migration_rate})')
    plt.yscale('log')
    plt.show()


def main():
    migration_rates = np.logspace(-6,0.5,9)  # replace with a list of your migration rates
    cutoff = 10**3 # replace with your chosen cutoff value

    probabilities = []
    errors = []

    for migration_rate in migration_rates:
        prob, err = calculate_probability_of_survival(migration_rate, cutoff)
        probabilities.append(prob)
        errors.append(err)

        # plot_total_infected_trajectories(migration_rate)  # Calling the function to plot the trajectories for each migration rate


    # Plotting the results
    plt.errorbar(migration_rates, probabilities, yerr=errors, fmt='o-', capsize=5)
    plt.xscale('log')
    plt.xlabel('Migration Rate')
    plt.ylabel('Probability of Survival')
    plt.title('Probability of Survival as a Function of Migration Rate')
    plt.show()

if __name__ == '__main__':
    main()
