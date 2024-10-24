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

migration_rates = 10**(np.linspace(-7, -0.5, 10))
N0 = 100
D = 1/100

# Read the CSV file
df = pd.read_csv("analysis_results.csv")

# Define colors for horizontal lines and ribbons
color_first = 'royalblue'
color_last = 'darkorange'

# Define extended migration rate range for shading
migration_rate_min = 1e-12
migration_rate_max = 10

# Calculate the standard error for a binomial distribution
n = 10000  # number of replicates
df['StandardError'] = np.sqrt(df['SurvivalProbability'] * (1 - df['SurvivalProbability']) / n)

# Extracting the first, last, and middle elements
first = df.iloc[0]
last = df.iloc[-1]
middle_df = df.iloc[1:10]

# New data
df_new = pd.read_csv("analysis_results_noback.csv")
n_new = 10000  # Assuming the same number of replicates for the new data
df_new['StandardError'] = np.sqrt(df_new['SurvivalProbability'] * (1 - df_new['SurvivalProbability']) / n_new)
middle_df_new = df_new.iloc[1:10]  # Adjust as per your new data structure

## TODO 
# 1. figsize = 6.5, 2.5 per row tall (this is in inches)
# 2. axis labels = 12
# 3. legends, insets = 10
# 4. subplot labels should be unbolded (a), (b), 

# Set the global font to Times New Roman
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 12  # Default font size
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 10
written_text_fontsize = 12

# Begin multi-panel figure setup
# fig, axs = plt.subplots(3, 2, figsize=(24, 16))  # Adjust overall figure size as needed
fig, axs = plt.subplots(3, 2, figsize=(6.5, 7.5))  # Adjust the figure size to balance the subplot shapes
labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

for ax, label in zip(axs.flat, labels):
    # Position: x=0, y=1 in axis coordinates (top left corner), transform=ax.transAxes
    # uses the axes coordinate system (0,0 is bottom left of the plot and 1,1 is top right)
    ax.text(0.02, 0.98, label, transform=ax.transAxes, fontweight='bold', va='top')

plt.subplots_adjust(wspace=0.3, hspace=0.3)  # Adjust spacing as needed

# Panel (f) plotting adjusted to fit into the subplot grid
ax = axs[0, 1]  # This selects the last panel position in a 3x2 grid

# Plotting code adjusted for the subplot
ax.errorbar(middle_df['MigrationRate'], middle_df['SurvivalProbability'], yerr=middle_df['StandardError'], fmt='-o', capsize=5, color='saddlebrown', label='Two-way Migration (1 \u2194 2)')
ax.errorbar(middle_df_new['MigrationRate'], middle_df_new['SurvivalProbability'], yerr=middle_df_new['StandardError'], fmt='-v', capsize=5, color='darkgreen', label='One-way Migration (1 \u2192 2)')


ax.axhline(y=first['SurvivalProbability'], color='royalblue', linestyle='--', label='Demixed (1 Only)')

ax.fill_between([migration_rate_min, migration_rate_max], first['SurvivalProbability'] - first['StandardError'], first['SurvivalProbability'] + first['StandardError'], color=color_first, alpha=0.3)

ax.axhline(y=last['SurvivalProbability'], color='darkorange', linestyle='--', label='Well-mixed (1 + 2)')
ax.fill_between([migration_rate_min, migration_rate_max], last['SurvivalProbability'] - last['StandardError'], last['SurvivalProbability'] + last['StandardError'], color=color_last, alpha=0.3)

ax.set_xscale('log')
ax.legend()
ax.set_xlim(1e-8, 1)
ax.set_ylim(bottom=0)
ax.set_xlabel(r'Migration Rate (Units: $\gamma$)')
ax.set_ylabel('Survival Probability')
ax.tick_params(axis='x', pad=15)
# ax.set_title('Survival Probability as a Function of Migration Rate')


# Panel (a) plotting adjusted for the subplot
ax = axs[0, 0]  # This selects the first panel position in a 3x2 grid

trajectories_df = pd.read_csv("trajectories_migration_rate_idx_1.csv")

# Adjusting code for panel (a)
for index, row in trajectories_df.iterrows():
    trajectory = ast.literal_eval(row['parent'])  # Convert string representation of list to actual list
    ax.plot(trajectory, linewidth=1.5)

start_time = 65
end_time = max([len(traj) for traj in trajectories_df['parent'].apply(ast.literal_eval)])

# Assuming the commented out lines are to be included, uncomment these as needed
# ax.axvspan(start_time, end_time, color='green', alpha=0.3)
ax.text(x=start_time + (end_time - start_time) / 2.3, y=1e6, s="Survival", rotation='horizontal', ha='center', va='center', color='green', fontsize=written_text_fontsize)

# ax.axhspan(0, 2, color='red', alpha=0.3)
ax.text(x=end_time * 0.74, y=5, s="Extinction", ha='center', va='center', color='red', fontsize=written_text_fontsize)

ax.set_yscale('log')
ax.set_xlabel(r'Time (Units: $\gamma^{-1}$)')
ax.set_ylabel('Total Infected Number')
ax.set_xlim(0, 80)
ax.set_ylim(bottom=1)
x_pos = np.max(ax.get_xlim()) * 0.95  # Adjust based on your axis limits
y_pos = np.min(ax.get_ylim()) * 3  # Pick a position within the visible range on a log scale

# # Plotting the dot
# ax.plot(x_pos, y_pos, 'o', color='blue', markersize=50, alpha=0.2, markeredgecolor='black')

# # Adding text in the middle of the dot
# ax.text(x_pos, y_pos, '1', color='black', verticalalignment='center', horizontalalignment='center')


# Panel (b) setup and plotting within the multi-panel figure
ax = axs[1, 0]  # Accessing the subplot for panel (b)

# Function for standard error calculation
def standard_error_of_proportion(p, n):
    return np.sqrt(p * (1 - p) / n)

# Load the data
csv_file_b = "special_case_antigenic_variance_deme1_migration_rate_idx_1.csv"  # Adjust the path as needed
data_b = pd.read_csv(csv_file_b)

# Filter and binning setup
xmin, xmax = 0.1, 0.3  # Adjust this range as needed
data_b = data_b[(data_b['AntigenicVariance'] >= xmin) & (data_b['AntigenicVariance'] <= xmax)]
num_bins = 30  # Number of bins
hist, bin_edges = np.histogram(data_b['AntigenicVariance'], bins=num_bins, range=(xmin, xmax), density=True)

# Calculate survival proportions and standard errors
survival_counts = []
non_survival_counts = []
standard_errors = []
for i in range(len(bin_edges)-1):
    bin_mask = (data_b['AntigenicVariance'] >= bin_edges[i]) & (data_b['AntigenicVariance'] < bin_edges[i+1])
    bin_data = data_b[bin_mask]
    n = len(bin_data)
    num_survived = bin_data['Survival'].sum()
    survival_counts.append(num_survived)
    non_survival_counts.append(n - num_survived)
    survival_proportion = num_survived / n if n > 0 else np.nan  # Avoid division by zero
    err = standard_error_of_proportion(survival_proportion, n) if n > 0 else np.nan
    standard_errors.append(err)

# Normalize the counts to get proportions
total_counts = np.array(survival_counts) + np.array(non_survival_counts)
survival_proportions = np.array(survival_counts) / total_counts
non_survival_proportions = np.array(non_survival_counts) / total_counts

# Plotting for panel (b)
widths = np.diff(bin_edges)
ax.bar(bin_edges[:-1], hist, width=widths, color='#87CEEB', edgecolor='none', align='edge')
ax.set_xlabel(r'Antigenic Diversity at Outbreak Max, $V(T_1)$')
ax.set_ylabel('Probability Density')
ax.set_xlim(xmin, xmax)
ax.tick_params(axis='x', pad=15)
# ax.tick_params(axis='y', pad=10)
# ax.set_title('Histogram of Antigenic Diversity & Probability of Survival')

# Define your color
curve_color = '#C71585'

# Second y-axis for probability of survival
ax2 = ax.twinx()
ax2.errorbar(bin_edges[:-1] + widths / 2, survival_proportions, yerr=standard_errors, color=curve_color, marker='o', linestyle='-', capsize=5)
ax2.set_ylabel('Probability of Survival', color=curve_color)  # Set y-axis label color

# Set the color of the tick labels
ax2.tick_params(axis='y', colors=curve_color)

# Optionally, set the color of the spine (the axis line)
ax2.spines['right'].set_color(curve_color)

# Mask to ignore NaN values for calculations
valid_mask = ~np.isnan(survival_proportions)

# Calculate the upper limit for y-axis considering only valid (non-NaN) survival_proportions
upper_ylim = np.max(survival_proportions[valid_mask] + np.array(standard_errors)[valid_mask])

# Set a bit higher ylim to ensure everything fits nicely
ax2.set_ylim(0, upper_ylim * 1.1)


from scipy.stats import poisson

# Setting up sample data for Poisson distribution for panel (c)
mu = 70  # Mean of distribution
x = np.arange(0, 200, 1)  # Discrete interval
pmf = poisson.pmf(x, mu)  # Poisson PMF for original and shifted data

# Accessing the subplot for panel (c)
ax = axs[2, 1]  # Position for panel (c) within the 3x2 grid

df = pd.read_csv("single_trajectory.csv")

# Plot total infected for each deme
ax.plot(df['times'], df['total_infected_1'], label='Deme 1', color='blue', linewidth=2)
ax.plot(df['times'], df['total_infected_2'], label='Deme 2', color='green', linewidth=2)

# Create a second y-axis for antigenic variance
ax2 = ax.twinx()
ax2.plot(df['times'], df['antigenic_variance_1'], label='Antigenic Variance Deme 1', color='blue', linestyle='-.', linewidth=2)
ax2.plot(df['times'], df['antigenic_variance_2'], label='Antigenic Variance Deme 2', color='green', linestyle='-.', linewidth=2)

max_idx_deme1 = df['total_infected_1'].idxmax()
max_idx_deme2 = df['total_infected_2'].idxmax()

time_max_infected_deme1 = df.loc[max_idx_deme1, 'times']
time_max_infected_deme2 = df.loc[max_idx_deme2, 'times']

# Vertical dashed lines
ax.axvline(time_max_infected_deme1, color='blue', linestyle='--', linewidth=2)
ax.axvline(time_max_infected_deme2, color='green', linestyle='--', linewidth=2)

# Annotate the first line
# ax.annotate('Max Infected Deme 1', xy=(time_max_infected_deme1, ax.get_ylim()[1]), xytext=(time_max_infected_deme1, ax.get_ylim()[1]*1.05),
#              arrowprops=dict(facecolor='blue', shrink=0.05), horizontalalignment='right')

# Draw a two-headed arrow between the lines
ax.annotate('', xy=(time_max_infected_deme1, ax.get_ylim()[1]*0.95), xytext=(time_max_infected_deme2, ax.get_ylim()[1]*0.95),
            arrowprops=dict(arrowstyle="<->", color='black'))

# Label the arrow as \Delta
mid_point = (time_max_infected_deme1 + time_max_infected_deme2) / 2
ax.text(mid_point, ax.get_ylim()[1]*0.86, r'$\Delta$', horizontalalignment='center', color='black')
ax.text(time_max_infected_deme1 - 0.9, ax.get_ylim()[1]*0.9, r'$T_1$', color='black', horizontalalignment='center')

ax2.plot(df.loc[max_idx_deme1, 'times'], df.loc[max_idx_deme1, 'antigenic_variance_1'], 'o', color='blue', markersize=15)
ax2.plot(df.loc[max_idx_deme2, 'times'], df.loc[max_idx_deme2, 'antigenic_variance_2'], 'o', color='green', markersize=15)

# Find the antigenic variances at the times of maximum infection
antigenic_variance_at_max1 = df.loc[max_idx_deme1, 'antigenic_variance_1']
antigenic_variance_at_max2 = df.loc[max_idx_deme2, 'antigenic_variance_2']

# Draw horizontal dashed lines from each dot to just before the secondary y-axis
# Note: ax2.get_xlim()[1] might need to be adjusted if it does not exactly match the secondary y-axis's end
ax2.plot([time_max_infected_deme1, ax2.get_xlim()[1]], [antigenic_variance_at_max1, antigenic_variance_at_max1], 'blue', linestyle='--', linewidth=2, dashes=(5, 5))
ax2.plot([time_max_infected_deme2, ax2.get_xlim()[1]], [antigenic_variance_at_max2, antigenic_variance_at_max2], 'green', linestyle='--', linewidth=2, dashes=(5, 5))

# Labeling
# ax.set_xlabel('Time')
ax.set_ylabel('Total Infected Number')
ax2.set_ylabel('Antigenic Diversity')
ax.set_ylim(bottom=0)
ax2.set_ylim(0, 1)

# Add legends
ax.legend(loc='center left')

# Labeling axes
ax.set_xlabel(r'Time (Units: $\gamma^{-1}$)')
ax.set_ylabel('Total Infected Number')
ax.set_xlim(0,22)

# Panel (e) setup for averages plotting, adapted for the subplot
ax = axs[1, 1]  # Accessing the subplot for panel (d)

cmap = plt.get_cmap('viridis')
colors = cmap(np.linspace(0, 1, len(migration_rates)))

avg_x_data = []
avg_y_data = []

for idx, migration_rate in enumerate(migration_rates):
    # Adjust file paths as necessary
    peak_time_path = f"peak_time_difference_migration_rate_idx_{idx + 2}.csv"
    peak_time_data = pd.read_csv(peak_time_path)['PeakTimeDifference'].values

    variance_diff_path = f"variance_difference_migration_rate_idx_{idx + 2}.csv"
    variance_diff_data = pd.read_csv(variance_diff_path)['VarianceDifference'].values

    avg_x_data.append(np.mean(2 * D * peak_time_data))
    avg_y_data.append(np.mean(variance_diff_data))

    migration_rate_label = f"{migration_rate:.1e}"
    ax.scatter(avg_x_data[-1], avg_y_data[-1], color=colors[idx], edgecolor='black', s=100, label=migration_rate_label)

# Plot theoretical curve (adjust as needed based on your theory or model)
xs = np.linspace(0, max(avg_x_data), 200) / 2 / D
ys = 2 * D * xs * (1 - np.exp(-20/2/N0 * (10 - xs)))
# ax.plot(2 * D * xs, ys)

# Plot y = x line
max_limit = max(max(avg_x_data), max(avg_y_data))
ax.plot([0, max_limit], [0, max_limit], 'k--', label='y = x')  # Black dashed line

ax.set_xlabel('2D * Average Peak Time Difference')
ax.set_ylabel('Average Variance Difference')
# ax.set_title('Averages of 2D * Peak Time Difference vs. Variance Difference')
# ax.legend(title="Migration Rate")
ax.grid(True, which="both", linestyle='--', linewidth=0.5)

# Panel (d) setup for average peak time difference plotting
ax = axs[2, 0]  # Accessing the subplot for panel (d), adjusted based on your correction

num_migration_rates = len(migration_rates)
average_peak_time_differences = []

for migration_rate_idx in range(2, num_migration_rates + 2):
    file_path = f"peak_time_difference_migration_rate_idx_{migration_rate_idx}.csv"
    peak_time_data = pd.read_csv(file_path)['PeakTimeDifference'].values
    
    average_peak_time = np.mean(peak_time_data)
    average_peak_time_differences.append(average_peak_time)

# Plotting directly within the subplot
ax.plot(migration_rates, average_peak_time_differences, marker='o', linestyle='-', color='blue')

ax.set_xlabel(r'Migration Rate (Units: $\gamma^{-1}$)')
ax.set_ylabel(r'Average Peak Time Difference, $\Delta$ (Units: $\gamma^{-1}$)')
ax.set_ylim(bottom=0)
# ax.set_title('Average Peak Time Difference vs. Migration Rate')
ax.set_xscale('log')
ax.tick_params(axis='x', pad=15)

# ax.set_yscale('log')
# ax.grid(True, which="both", linestyle='--', linewidth=0.5)

from scipy.integrate import trapz

# Given parameters
N0 = 100
sig = np.sqrt(20)
F = 1.5
D = 0.01
T1 = 10

# CDF Delta function
def cdf_delta(t, k):
    return np.exp(-2 * N0 * k / sig**2 * (np.exp(F * t) - 1))

# Arrays for t and k
ts = np.linspace(0, 15, int(1e4))
ks = np.logspace(-8, -0.5, 20)

# Calculate avg_deltas
avg_deltas = np.array([trapz(cdf_delta(ts, k), ts) for k in ks])

# Assuming migration_rates are equivalent to ks here
migration_rates = ks

ax.plot(migration_rates, avg_deltas, linestyle='--', color='red', label='Theoretical Estimate')
# ax.plot(migration_rates, 1/F * np.log(20 / 2 / N0 / migration_rates), linestyle='--', color='red')  # Theoretical curve, adjust color/style as needed








# Show the entire multi-panel figure
# plt.tight_layout(pad=3.0, h_pad=3.0, w_pad=2.0)  # Adjust padding and spacing as needed
fig.set_constrained_layout(True)

plt.savefig('twodeme_multipanel_coded.pdf', format='pdf', dpi=300)
plt.savefig('twodeme_multipanel_coded.png', format='png')
plt.savefig('twodeme_multipanel_coded.svg', format='svg', dpi=300)

# plt.show()