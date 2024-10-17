import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Ellipse
import numpy as np
import ast  # Importing the ast module

# Read the CSV file
df = pd.read_csv("three_deme_one_edge.csv")

# Calculate the standard error for a binomial distribution
n = 10000  # number of replicates
df['StandardError'] = np.sqrt(df['SurvivalProbability'] * (1 - df['SurvivalProbability']) / n)

# Extracting the first, last, and middle elements
first = df.iloc[0]
# last = df.iloc[-1]
middle_df = df.iloc[0:10]

# Set the global font to Times New Roman
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 14  # Default font size
rcParams['axes.labelsize'] = 16
rcParams['axes.titlesize'] = 18
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
rcParams['legend.fontsize'] = 16

# Plotting
plt.figure(figsize=(10, 6))

# Plot middle elements with error bars
plt.errorbar(middle_df['MigrationRate'], middle_df['SurvivalProbability'], yerr=middle_df['StandardError'], fmt='o', capsize=5, label='Simulation')

# Plot first and last as horizontal lines
plt.axhline(y=first['SurvivalProbability'], color='r', linestyle='--', label='Demixed')
# plt.axhline(y=last['SurvivalProbability'], color='g', linestyle='--', label='Well-mixed')

plt.xscale('log')
plt.legend()
plt.xlabel('Migration Rate')
plt.ylabel('Survival Probability')
plt.title('Survival Probability as a Function of Migration Rate')

plt.savefig('Survival_Probability_three_deme_one_edge.pdf', format='pdf', dpi=300)
plt.savefig('Survival_Probability_three_deme_one_edge.png', format='png')

