import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

# Set the global font to Times New Roman
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 12  # Default font size
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 9

# Load the CSV file
df = pd.read_csv("survival_probabilities.csv")

# Define migration rates and host population sizes
migration_rates = np.logspace(-6.0, -2, 30)
host_population_sizes = np.logspace(4.0, 7.0, 7)

# Find the optimal migration rate for each host population size
opt_migration = []
for k in range(len(host_population_sizes)):
    subset = df[df['HostPerDemeIdx'] == k]
    if len(subset) > 0:
        y_data = subset['SurvivalProbability'].values
        max_idx = np.argmax(y_data)
        opt_migration.append((host_population_sizes[k], migration_rates[max_idx]))

# Extract host sizes and optimal migration rates for plotting
host_sizes, optimal_migration_rates = zip(*opt_migration)

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(host_sizes, optimal_migration_rates, color='b', label=None)

# Add dashed line showing migration proportional to host size ^ -1/2
host_sizes_line = np.logspace(4.0, 7.0, 100)
migration_line = host_sizes_line ** (-0.5)
migration_line_2 = host_sizes_line ** (-1)
plt.plot(host_sizes_line * host_sizes[0] / host_sizes_line[0], migration_line * optimal_migration_rates[0] / migration_line[0], linestyle='--', color='r', label=r'$k \propto N_h^{-1/2}$')
plt.plot(host_sizes_line * host_sizes[0] / host_sizes_line[0], migration_line_2 * optimal_migration_rates[0] / migration_line_2[0], linestyle='--', color='k', label=r'$k \propto N_h^{-1}$')

# Set the x-axis and y-axis to log scale
plt.xscale('log')
plt.yscale('log')

# Add labels and title
plt.xlabel(r'host population size $N_h$')
plt.ylabel(r'optimal migration rate, $k^*/\gamma$')
# plt.title('Optimal Migration Rate vs Host Population Size')
plt.legend(loc='upper right', framealpha=0.3)

plt.tight_layout()

# Display the plot
# plt.show()

# Optionally, save the plot to a file
plt.savefig("peak_migration_rate_vs_host_population_size.png")
plt.savefig("peak_migration_rate_vs_host_population_size.pdf")
plt.savefig("peak_migration_rate_vs_host_population_size.svg")
