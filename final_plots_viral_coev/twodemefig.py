import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Circle
import ast  # Importing the ast module

# Read the CSV file
df = pd.read_csv("analysis_results.csv")

# Calculate the standard error for a binomial distribution
n = 10000  # number of replicates
df['StandardError'] = np.sqrt(df['SurvivalProbability'] * (1 - df['SurvivalProbability']) / n)

# Extracting the first, last, and middle elements
first = df.iloc[0]
last = df.iloc[-1]
middle_df = df.iloc[1:10]

# Set the global font to Times New Roman
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 20  # Default font size
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 18

# Plotting
plt.figure(figsize=(12, 8))

# Plot middle elements with error bars and lines
plt.errorbar(middle_df['MigrationRate'], middle_df['SurvivalProbability'], yerr=middle_df['StandardError'], fmt='-o', capsize=5, color='saddlebrown', label='Simulation')

# Define colors for horizontal lines and ribbons
color_first = 'royalblue'
color_last = 'darkorange'

# Define extended migration rate range for shading
migration_rate_min = 1e-12
migration_rate_max = 10

# Plot first and last as horizontal lines with extended ribbons
plt.axhline(y=first['SurvivalProbability'], color=color_first, linestyle='--', label='Demixed')
plt.fill_between([migration_rate_min, migration_rate_max], first['SurvivalProbability'] - first['StandardError'], first['SurvivalProbability'] + first['StandardError'], color=color_first, alpha=0.3)

plt.axhline(y=last['SurvivalProbability'], color=color_last, linestyle='--', label='Well-mixed')
plt.fill_between([migration_rate_min, migration_rate_max], last['SurvivalProbability'] - last['StandardError'], last['SurvivalProbability'] + last['StandardError'], color=color_last, alpha=0.3)

plt.xscale('log')
plt.legend()
plt.xlim(1e-8, 1)
plt.ylim(bottom=0)
plt.xlabel('Migration Rate')
plt.ylabel('Survival Probability')
plt.title('Survival Probability as a Function of Migration Rate')

plt.savefig('Survival_Probability.pdf', format='pdf', dpi=300)
plt.savefig('Survival_Probability.png', format='png')

# Read the trajectories CSV file
# Replace 'trajectories.csv' with your actual file name
trajectories_df = pd.read_csv("trajectories_migration_rate_idx_1.csv")

# Plotting Trajectories
plt.figure(figsize=(12, 8))

# Convert string representations of lists to actual lists and plot each trajectory
for index, row in trajectories_df.iterrows():
    trajectory = ast.literal_eval(row['parent'])  # Safely evaluate the string as a list
    plt.plot(trajectory, linewidth=1.5)

# Shade all time points from 75 onward
start_time = 65
end_time = max([len(traj) for traj in trajectories_df['parent'].apply(ast.literal_eval)])
# plt.axvspan(start_time, end_time, color='green', alpha=0.3)
# Add vertical text "Survival" within the shaded area
plt.text(x=start_time + (end_time - start_time) / 2, y=1e6, s="Survival", rotation='horizontal', ha='center', va='center', color='green', fontsize=20)

# Add a red horizontal rectangle labeled "Extinction"
# plt.axhspan(0, 2, color='red', alpha=0.3)  # y-range slightly around 1
plt.text(x=end_time * 0.67, y=1, s="Extinction", ha='center', va='center', color='red', fontsize=20)

plt.yscale('log')
plt.xlabel('Time')
plt.ylabel('Total Infected Number')
plt.xlim(0, 80)
plt.title('Total Infected Number vs Time for 200 Trajectories')

plt.savefig('Infected_Trajectories.pdf', format='pdf', dpi=300)
plt.savefig('Infected_Trajectories.png', format='png')

# Function to calculate standard error of proportion
def standard_error_of_proportion(p, n):
    return np.sqrt(p * (1 - p) / n)

# Load the data
csv_file = "special_case_antigenic_variance_deme1_migration_rate_idx_1.csv"  # Adjust the path as needed
data = pd.read_csv(csv_file)

# Define the range for x-axis and binning
xmin, xmax = 0.1, 0.3  # Adjust this range as needed
data = data[(data['AntigenicVariance'] >= xmin) & (data['AntigenicVariance'] <= xmax)]

# Calculate the histogram within the specified range
num_bins = 15  # Number of bins
hist, bin_edges = np.histogram(data['AntigenicVariance'], bins=num_bins, range=(xmin, xmax), density=False)

# Calculate the proportion of surviving and non-surviving trajectories in each bin
survival_counts = []
non_survival_counts = []
standard_errors = []

for i in range(len(bin_edges)-1):
    bin_mask = (data['AntigenicVariance'] >= bin_edges[i]) & (data['AntigenicVariance'] < bin_edges[i+1])
    bin_data = data[bin_mask]
    n = len(bin_data)
    num_survived = bin_data['Survival'].sum()
    survival_counts.append(num_survived)
    non_survival_counts.append(n - num_survived)
    survival_proportion = num_survived / n
    standard_errors.append(standard_error_of_proportion(survival_proportion, n))

# Normalize the counts to get proportions
total_counts = np.array(survival_counts) + np.array(non_survival_counts) 
survival_proportions = np.array(survival_counts) / total_counts
non_survival_proportions = np.array(non_survival_counts) / total_counts

# Plotting
fig, ax1 = plt.subplots(figsize=(12, 8))

# Warm, colorblind-safe colors
colors_survived = '#87CEEB'  # Sky blue
colors_not_survived = '#FFD700'  # Yellow

# Stacked bar plot
widths = np.diff(bin_edges)
ax1.bar(bin_edges[:-1], hist, width=widths, color=colors_survived, edgecolor='gray', align='edge')
# ax1.bar(bin_edges[:-1], survival_proportions * hist, width=widths, color=colors_survived, edgecolor='gray', align='edge', label='Survived')
# ax1.bar(bin_edges[:-1], non_survival_proportions * hist, bottom=survival_proportions * hist, width=widths, color=colors_not_survived, edgecolor='gray', align='edge', label='Not Survived')
ax1.set_xlabel('Antigenic Diversity at Outbreak Max')
ax1.set_ylabel('Counts')
ax1.set_xlim(xmin, xmax)
ax1.set_title('Histogram of Antigenic Diversity & Probability of Survival')

# Add second y-axis for the probability of survival
ax2 = ax1.twinx()
ax2.errorbar(bin_edges[:-1] + widths / 2, survival_proportions, yerr=standard_errors, color='#C71585', marker='o', linestyle='-', label='Probability of Survival', capsize=5)  # Dark Fuchsia
ax2.set_ylabel('Probability of Survival')
ax2.legend(loc='upper right')

plt.savefig('singledeme_variance_histogram.pdf', format='pdf', dpi=300)
plt.savefig('singledeme_variance_histogram.png', format='png')

# Plotting
plt.figure(figsize=(12, 8))

num_bins = 15

# Overall Histogram (Stairstep)
hist, bin_edges, _ = plt.hist(data['AntigenicVariance'], bins=num_bins, range=(xmin, xmax), density=True, histtype='step', color='black', label='Overall')

# Histogram Conditioned on Survival
survived_data = data[data['Survival'] == 1]
plt.hist(survived_data['AntigenicVariance'], bins=bin_edges, density=True, alpha=0.5, color='#87CEEB', label='Survived')

# Histogram Conditioned on Non-Survival
# non_survived_data = data[data['Survival'] == 0]
# plt.hist(non_survived_data['AntigenicVariance'], bins=bin_edges, density=True, alpha=0.5, color='#FFD700', label='Not Survived')

plt.xlabel('Antigenic Diversity at Outbreak Max')
plt.ylabel('Density')
plt.title('Histogram of Antigenic Diversity')
plt.xlim(xmin, xmax)
plt.legend()
plt.savefig('survived_histogram.pdf', format='pdf', dpi=300)
plt.savefig('survived_histogram.png', format='png')

# Existing data
df = pd.read_csv("analysis_results.csv")
n = 10000  # number of replicates for the existing data
df['StandardError'] = np.sqrt(df['SurvivalProbability'] * (1 - df['SurvivalProbability']) / n)
middle_df = df.iloc[1:10]

# New data
df_new = pd.read_csv("analysis_results_noback.csv")
n_new = 10000  # Assuming the same number of replicates for the new data
df_new['StandardError'] = np.sqrt(df_new['SurvivalProbability'] * (1 - df_new['SurvivalProbability']) / n_new)
middle_df_new = df_new.iloc[1:10]  # Adjust as per your new data structure

# Plot settings
plt.figure(figsize=(12, 8))

# Plot existing data points with the updated label
plt.errorbar(middle_df['MigrationRate'], middle_df['SurvivalProbability'], yerr=middle_df['StandardError'], fmt='-o', capsize=5, color='saddlebrown', label='Two-way Migration (1 \u2194 2)')

# Plot new data points with the updated label
plt.errorbar(middle_df_new['MigrationRate'], middle_df_new['SurvivalProbability'], yerr=middle_df_new['StandardError'], fmt='-v', capsize=5, color='darkgreen', label='One-way Migration (1 \u2192 2)')

# Customize plot
plt.xscale('log')
plt.xlim(1e-8, 1)
plt.ylim(bottom=0)
plt.xlabel('Migration Rate')
plt.ylabel('Survival Probability')
plt.title('Survival Probability as a Function of Migration Rate')

# Custom legend for horizontal lines (if they exist in your plot)
# Assuming first and last are defined as in your setup for horizontal line plotting
plt.axhline(y=first['SurvivalProbability'], color='royalblue', linestyle='--', label='Demixed (1 Only)')
plt.axhline(y=last['SurvivalProbability'], color='darkorange', linestyle='--', label='Well-mixed (1 + 2)')

# Retrieve handles and labels
handles, labels = plt.gca().get_legend_handles_labels()

# Suppose you want the order: Unidirectional, Bidirectional, Single entity, Combined entity
# Reorder handles and labels according to desired order
new_order = [2, 3, 1, 0]  # New order as indices of the original handles/labels
ordered_handles = [handles[idx] for idx in new_order]
ordered_labels = [labels[idx] for idx in new_order]

# Create legend with reordered handles and labels
plt.legend(ordered_handles, ordered_labels)

plt.savefig('survival_prob_w_noback.pdf', format='pdf', dpi=300)
plt.savefig('survival_prob_w_noback.png', format='png')


plt.show()