import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
import numpy as np
import ast
from scipy.stats import norm
from scipy.stats import linregress

# Set the global font to Times New Roman
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 12  # Default font size
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 9
written_text_fontsize = 12

# Load the CSV files
df_deme1 = pd.read_csv("probability_deme1_first.csv")
df_deme2 = pd.read_csv("probability_deme2_first.csv")
df_both = pd.read_csv("probability_both_same_time.csv")

# Extract data
migration_rates = df_deme1['MigrationRate'].values
prob_deme1_first = df_deme1['ProbabilityDeme1First'].values
prob_deme2_first = df_deme2['ProbabilityDeme2First'].values
prob_both_same_time = df_both['ProbabilityBothSameTime'].values

# Create the stacked plot
plt.figure(figsize=(6.5/2, 2.75))
plt.fill_between(migration_rates, 0, prob_both_same_time, color='green', alpha=0.5, label='both at same time')
plt.fill_between(migration_rates, prob_both_same_time, prob_both_same_time + prob_deme2_first, color='red', alpha=0.5, label='deme 2 first')
plt.fill_between(migration_rates, prob_both_same_time + prob_deme2_first, prob_both_same_time + prob_deme2_first + prob_deme1_first, color='blue', alpha=0.5, label='deme 1 first')

# Set the x-axis to log scale
plt.xscale('log')

# Add labels and title
plt.xlabel(r'migration rate, $k/\gamma$')
plt.ylabel('prob. of first escape strain')
# plt.title('Probability of Outbreak Occurrence by Migration Rate')
plt.legend(loc='upper left', framealpha=0.4)
plt.xlim(migration_rates[0], migration_rates[-1])
plt.ylim(0,1)

plt.tight_layout()

# Display the plot
# plt.show()

# Optionally, save the plot to a file
plt.savefig("cumulative_outbreak_probabilities_migration_rate.png")
plt.savefig("cumulative_outbreak_probabilities_migration_rate.pdf")
plt.savefig("cumulative_outbreak_probabilities_migration_rate.svg")


