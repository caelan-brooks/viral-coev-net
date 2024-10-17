import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from matplotlib import rcParams
import matplotlib.ticker as ticker

# Set the global font to Times New Roman
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 12  # Default font size
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 12

# Load the data
csv_file = "survival_probabilities.csv"  # Adjust the file name/path if necessary
data = pd.read_csv(csv_file)

# Define the migration and mutation rates as in the Julia script
migration_rates = np.power(10.0, np.linspace(-7, -0.5, 10))
mutation_rates = np.linspace(0.001, 0.02, 10)

# Pivot the DataFrame for the heatmap
pivot_table = data.pivot("MigrationRateIdx", "MutationRateIdx", "SurvivalProbability")

# Create the heatmap
plt.figure(figsize=(10, 8))
ax = plt.gca()  # Get current axis
cax = ax.imshow(pivot_table.values, interpolation='nearest', cmap='viridis', aspect='auto')

# Add color bar
cbar = plt.colorbar(cax)
cbar.set_label('Survival Probability')

# Formatting tick labels in scientific notation with one decimal place
ax.set_xticks(np.arange(len(mutation_rates)))
ax.set_xticklabels([f'{rate:.1e}' for rate in mutation_rates])
ax.set_yticks(np.arange(len(migration_rates)))
ax.set_yticklabels([f'{rate:.1e}' for rate in migration_rates])

# Adding labels
plt.title("Heatmap of Survival Probability vs Mutation and Migration Rates")
plt.xlabel("Mutation Rate")
plt.ylabel("Migration Rate")

plt.savefig('mutation_migration_heatmap.pdf', format='pdf', dpi=300)
plt.savefig('mutation_migration_heatmap.png', format='png')

# Function to calculate standard error
def standard_error(p, n=10000):
    return np.sqrt(p * (1 - p) / n)

# Plotting survival probability vs migration rate for different mutation rates
plt.figure(figsize=(12, 8))
for mutation_rate_idx in range(1, len(mutation_rates) + 1):
    subset = data[data["MutationRateIdx"] == mutation_rate_idx]
    errors = subset["SurvivalProbability"].apply(lambda p: standard_error(p))
    plt.errorbar(migration_rates, subset["SurvivalProbability"], yerr=errors, label=f'Mutation Rate {mutation_rates[mutation_rate_idx - 1]:.1e}', fmt='-o', markersize=4, capsize=5, elinewidth=1, errorevery=1)
plt.xlabel("Migration Rate")
plt.ylabel("Survival Probability")
plt.title("Survival Probability vs Migration Rate for Different Mutation Rates")
plt.xscale("log")
plt.legend(title="Mutation Rates", loc='upper left')
plt.grid(True)
plt.savefig('mutation_migration_curves_with_error.pdf', format='pdf', dpi=300)
plt.savefig('mutation_migration_curves_with_error.png', format='png')

# Plotting survival probability vs mutation rate for different migration rates
plt.figure(figsize=(12, 8))
for migration_rate_idx in range(1, len(migration_rates) + 1):
    subset = data[data["MigrationRateIdx"] == migration_rate_idx]
    errors = subset["SurvivalProbability"].apply(lambda p: standard_error(p))
    plt.errorbar(mutation_rates, subset["SurvivalProbability"], yerr=errors, label=f'{migration_rates[migration_rate_idx - 1]:.1e}', fmt='-o', markersize=4, capsize=5, elinewidth=1, errorevery=1)



plt.xlabel(r"mutation rate, $D$")
plt.ylabel("escpape probability")
# plt.title("Survival Probability vs Mutation Rate for Different Migration Rates")
plt.xscale("linear")
plt.legend(title=r"migration rates, $k/\gamma$", loc='upper left')
plt.grid(True)
# Using StrMethodFormatter for x-axis labels
plt.gca().xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))

plt.savefig('migration_mutation_curves_with_error.pdf', format='pdf', dpi=300)
plt.savefig('migration_mutation_curves_with_error.png', format='png')
plt.show()