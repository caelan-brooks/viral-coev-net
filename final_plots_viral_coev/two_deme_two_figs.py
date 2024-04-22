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
from scipy.stats import norm
from scipy.stats import linregress

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
middle_df = df.iloc[1:11]

# New data
df_new = pd.read_csv("analysis_results_noback.csv")
n_new = 10000  # Assuming the same number of replicates for the new data
df_new['StandardError'] = np.sqrt(df_new['SurvivalProbability'] * (1 - df_new['SurvivalProbability']) / n_new)
middle_df_new = df_new.iloc[1:11]  # Adjust as per your new data structure

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
rcParams['legend.fontsize'] = 9
written_text_fontsize = 12

## First figure will have the trajectories and the survival probabilities 
fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.5))  # Adjust the figure size to balance the subplot shapes
labels = ['(a)', '(b)']

for ax, label in zip(axs.flat, labels):
    # Position: x=0, y=1 in axis coordinates (top left corner), transform=ax.transAxes
    # uses the axes coordinate system (0,0 is bottom left of the plot and 1,1 is top right)
    ax.text(-0.01, 1.05, label, transform=ax.transAxes, va='bottom')

ax = axs[0]
trajectories_df = pd.read_csv("trajectories_migration_rate_idx_1.csv")

# Adjusting code for panel (a)
for index, row in trajectories_df.iterrows():
    trajectory = ast.literal_eval(row['parent'])  # Convert string representation of list to actual list
    ax.plot(trajectory, linewidth=1.5)

start_time = 65
end_time = max([len(traj) for traj in trajectories_df['parent'].apply(ast.literal_eval)])

# Assuming the commented out lines are to be included, uncomment these as needed
# ax.axvspan(start_time, end_time, color='green', alpha=0.3)
ax.text(x=start_time + (end_time - start_time) / 3.0, y=1e6, s="survival", rotation='horizontal', ha='center', va='center', color='green', fontsize=written_text_fontsize)

# ax.axhspan(0, 2, color='red', alpha=0.3)
ax.text(x=end_time * 0.73, y=3, s="extinction", ha='center', va='center', color='red', fontsize=written_text_fontsize)

ax.set_yscale('log')
ax.set_xlabel(r'time (units: $\gamma^{-1}$)')
ax.set_ylabel('total infected number')
ax.set_xlim(0, 80)
ax.set_ylim(bottom=1)
x_pos = np.max(ax.get_xlim()) * 0.95  # Adjust based on your axis limits
y_pos = np.min(ax.get_ylim()) * 3  # Pick a position within the visible range on a log scale


ax = axs[1]  # This selects the last panel position in a 3x2 grid

# Plotting code adjusted for the subplot
ax.errorbar(middle_df['MigrationRate'], middle_df['SurvivalProbability'], yerr=middle_df['StandardError'], fmt='-o', capsize=5, color='saddlebrown', label='Two-way Migration (1 \u2194 2)')
ax.errorbar(middle_df_new['MigrationRate'], middle_df_new['SurvivalProbability'], yerr=middle_df_new['StandardError'], fmt='-v', capsize=5, color='darkgreen', label='One-way Migration (1 \u2192 2)')


ax.axhline(y=first['SurvivalProbability'], color='royalblue', linestyle='--')

ax.fill_between([migration_rate_min, migration_rate_max], first['SurvivalProbability'] - first['StandardError'], first['SurvivalProbability'] + first['StandardError'], color=color_first, alpha=0.3)

ax.axhline(y=last['SurvivalProbability'], color='darkorange', linestyle='--')
ax.fill_between([migration_rate_min, migration_rate_max], last['SurvivalProbability'] - last['StandardError'], last['SurvivalProbability'] + last['StandardError'], color=color_last, alpha=0.3)

ax.set_xscale('log')
ax.legend(frameon=False, loc=(0.1, 0.01))
ax.set_xlim(1e-8, 1)
ax.set_ylim(bottom=0)
ax.set_xlabel(r'migration rate (units: $\gamma$)')
ax.set_ylabel('survival probability')
# ax.tick_params(axis='x', pad=15)
# ax.set_title('Survival Probability as a Function of Migration Rate')

#######################
plt.tight_layout()
plt.savefig('twodeme_trajectories_and_probs.pdf', format='pdf', dpi=300)
print('Done')

########################

fig, axs = plt.subplots(2, 2, figsize=(6.5, 5.0))  # Adjust the figure size to balance the subplot shapes
labels = ['(a)', '(b)', '(c)', '(d)']

for ax, label in zip(axs.flat, labels):
    # Position: x=0, y=1 in axis coordinates (top left corner), transform=ax.transAxes
    # uses the axes coordinate system (0,0 is bottom left of the plot and 1,1 is top right)
    ax.text(-0.02, 1.05, label, transform=ax.transAxes, va='bottom')


ax = axs[0, 0]

# Parameters for the Gaussians
mean = 0
sigma1 = np.sqrt(1)  # Standard deviation for deme 1
sigma2 = np.sqrt(0.2)  # Standard deviation for deme 2

# Function to calculate FWHM
def fwhm(sigma):
    return 2 * np.sqrt(2 * np.log(2)) * sigma

# Generate x values
x = np.linspace(-5, 5, 1000)

# Gaussian distributions
y1 = norm.pdf(x, mean, sigma1)
y2 = norm.pdf(x, mean, sigma2)

# Scale for visibility if needed, or adjust plot limits or aspect
scaling_factor = 7
y1_scaled = y1 * scaling_factor
y2_scaled = y2 

# Plotting on the specified axis
ax.plot(x, y1_scaled, label='deme 1', linewidth=2, color='blue')
ax.plot(x, y2_scaled, label='deme 2', linewidth=2, color='green')

# Calculate FWHM for each distribution
fwhm1 = fwhm(sigma1)
fwhm2 = fwhm(sigma2)

# Height for arrows (half-maximum of the unscaled distributions)
height1 = max(y1) * scaling_factor / 2
height2 = max(y2) / 2

# Plot two-way arrows for FWHM
ax.annotate('', xy=(mean - fwhm1/2, height1), xytext=(mean + fwhm1/2, height1),
            arrowprops=dict(arrowstyle='<->', color='blue'))
ax.annotate('', xy=(mean - fwhm2/1.5, height2), xytext=(mean + fwhm2/1.5, height2),
            arrowprops=dict(arrowstyle='<->', color='green'))

# Text annotations for FWHM
ax.text(mean, height1 * 0.95, '$V_1$', ha='center', va='top', color='blue')
ax.text(mean, height2 * 0.80, r'$V_2$', ha='center', va='top', color='green')

# Additional plot settings
ax.legend(frameon=False, loc=(0.58, 0.7))
ax.set_xlabel(r'antigenic coordinate, $x$')
ax.set_ylabel(r'infected number density $n_i(x,t)$')
ax.set_yticks([])
ax.set_xticks([])

ax = axs[1, 0]
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
ax.set_xlabel(r'antigenic diversity at outbreak max, $V(T_1)$')
ax.set_ylabel('probability density')
ax.set_xlim(xmin, xmax)
# ax.tick_params(axis='x', pad=15)
# ax.tick_params(axis='y', pad=10)
# ax.set_title('Histogram of Antigenic Diversity & Probability of Survival')

# Define your color
curve_color = '#C71585'

# Second y-axis for probability of survival
ax2 = ax.twinx()
ax2.errorbar(bin_edges[:-1] + widths / 2, survival_proportions, yerr=standard_errors, color=curve_color, marker='o', linestyle='-', capsize=2, markersize=3)
ax2.set_ylabel('probability of survival', color=curve_color)  # Set y-axis label color


# Set the color of the tick labels
ax2.tick_params(axis='y', colors=curve_color)

# Optionally, set the color of the spine (the axis line)
ax2.spines['right'].set_color(curve_color)

# Mask to ignore NaN values for calculations
valid_mask = ~np.isnan(survival_proportions)

# Calculate the upper limit for y-axis considering only valid (non-NaN) survival_proportions
upper_ylim = np.max(survival_proportions[valid_mask] + np.array(standard_errors)[valid_mask])

# Set a bit higher ylim to ensure everything fits nicely
ax2.set_ylim(0, upper_ylim * 1.3)

# Calculate bin centers
bin_centers = bin_edges[:-1] + widths / 2

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(bin_centers[valid_mask], survival_proportions[valid_mask])
print(slope, intercept)
# Calculate R^2 and Pearson correlation
r_squared = r_value**2
pearson_corr = r_value
x_fit = np.linspace(bin_centers.min(), bin_centers.max(), 100)
y_fit = slope * x_fit + intercept
ax2.plot(x_fit, y_fit, 'k--', label=f'Fit: $R^2$={r_squared:.2f}, \n Pearson={pearson_corr:.2f}')
ax2.legend(frameon=False, loc=(0.40, 0.75))
#######
ax = axs[1,1]
cmap = plt.get_cmap('inferno')
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

ax.set_xlabel('2D * average peak time difference')
ax.set_ylabel('average diversity difference')
# ax.set_title('Averages of 2D * Peak Time Difference vs. Variance Difference')
# ax.legend(title="Migration Rate")
ax.grid(True, which="both", linestyle='--', linewidth=0.5)

ax = axs[0, 1]
df = pd.read_csv("single_trajectory.csv")

# Plot total infected for each deme
ax.plot(df['times'], df['total_infected_1'], color='blue', linewidth=2)
ax.plot(df['times'], df['total_infected_2'], color='green', linewidth=2)

# Create a second y-axis for antigenic variance
ax2 = ax.twinx()
ax2.plot(df['times'], df['antigenic_variance_1'], color='blue', linestyle='-.', linewidth=2)
ax2.plot(df['times'], df['antigenic_variance_2'], color='green', linestyle='-.', linewidth=2)

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
ax.annotate('', xy=(time_max_infected_deme1, ax.get_ylim()[1]*2), xytext=(time_max_infected_deme2, ax.get_ylim()[1]*2),
            arrowprops=dict(arrowstyle="<->", color='black'))

# Label the arrow as \Delta
mid_point = (time_max_infected_deme1 + time_max_infected_deme2) / 2
ax.text(mid_point, ax.get_ylim()[1]*1/2, r'$\Delta$', horizontalalignment='center', color='black')
ax.text(time_max_infected_deme1-0.4, 0.2 , r'$T_1$', color='black', horizontalalignment='center')

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
ax.set_ylabel('infected number')
ax2.set_ylabel('antigenic diversity')
ax.set_yscale('log')
ax.set_ylim(bottom=1)
ax2.set_ylim(0, 1)

# Add legends
ax.legend(loc='center left', frameon=False)

# Labeling axes
ax.set_xlabel(r'time (units: $\gamma^{-1}$)')
ax.set_ylabel('total infected number')
ax.set_xlim(0,22)
#######################
plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0) 
plt.savefig('twodeme_antigenic_diversity.pdf', format='pdf', dpi=300)
print('Done again')
#######################

## First figure will have the trajectories and the survival probabilities 
fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.5))  # Adjust the figure size to balance the subplot shapes
labels = ['(a)', '(b)']

for ax, label in zip(axs.flat, labels):
    # Position: x=0, y=1 in axis coordinates (top left corner), transform=ax.transAxes
    # uses the axes coordinate system (0,0 is bottom left of the plot and 1,1 is top right)
    ax.text(-0.01, 1.05, label, transform=ax.transAxes, va='bottom')

ax = axs[0]
trajectories_df = pd.read_csv("trajectories_migration_rate_idx_1.csv")

# Adjusting code for panel (a)
for index, row in trajectories_df.iterrows():
    trajectory = ast.literal_eval(row['parent'])  # Convert string representation of list to actual list
    ax.plot(trajectory, linewidth=1.5)

start_time = 65
end_time = max([len(traj) for traj in trajectories_df['parent'].apply(ast.literal_eval)])

# Assuming the commented out lines are to be included, uncomment these as needed
# ax.axvspan(start_time, end_time, color='green', alpha=0.3)
ax.text(x=start_time + (end_time - start_time) / 3.0, y=1e6, s="survival", rotation='horizontal', ha='center', va='center', color='green', fontsize=written_text_fontsize)

# ax.axhspan(0, 2, color='red', alpha=0.3)
ax.text(x=end_time * 0.73, y=3, s="extinction", ha='center', va='center', color='red', fontsize=written_text_fontsize)

ax.set_yscale('log')
ax.set_xlabel(r'time (units: $\gamma^{-1}$)')
ax.set_ylabel('total infected number')
ax.set_xlim(0, 80)
ax.set_ylim(bottom=1)
x_pos = np.max(ax.get_xlim()) * 0.95  # Adjust based on your axis limits
y_pos = np.min(ax.get_ylim()) * 3  # Pick a position within the visible range on a log scale


ax = axs[1]  # This selects the last panel position in a 3x2 grid

# Plotting code adjusted for the subplot
ax.errorbar(middle_df['MigrationRate'], middle_df['SurvivalProbability'], yerr=middle_df['StandardError'], fmt='-o', capsize=5, color='saddlebrown', label='Two-way Migration (1 \u2194 2)')
ax.errorbar(middle_df_new['MigrationRate'], middle_df_new['SurvivalProbability'], yerr=middle_df_new['StandardError'], fmt='-v', capsize=5, color='darkgreen', label='One-way Migration (1 \u2192 2)')


ax.axhline(y=first['SurvivalProbability'], color='royalblue', linestyle='--')

ax.fill_between([migration_rate_min, migration_rate_max], first['SurvivalProbability'] - first['StandardError'], first['SurvivalProbability'] + first['StandardError'], color=color_first, alpha=0.3)

ax.axhline(y=last['SurvivalProbability'], color='darkorange', linestyle='--')
ax.fill_between([migration_rate_min, migration_rate_max], last['SurvivalProbability'] - last['StandardError'], last['SurvivalProbability'] + last['StandardError'], color=color_last, alpha=0.3)

def load_and_compute_averages(num_migration_rates):
    avg_variances = []
    for idx in range(2, num_migration_rates + 2):
        deme2_data = pd.read_csv(f"antigenic_variance_deme2_migration_rate_idx_{idx}.csv")['AntigenicVariance'].values
        avg_variance = np.mean(deme2_data) * len(deme2_data) / 10001 + (-intercept/slope) * (1 - len(deme2_data)/10001) # adjust for the fact that at low migration not all trajectories result in secondary outbreaks
        avg_variances.append(avg_variance)
    return np.array(avg_variances)

# Calculate p2 and the modified survival probability
def compute_probabilities(avg_variances, p1, slope, intercept):
    p2_values = slope * avg_variances + intercept
    final_probabilities = 1 - (1 - p1) * (1 - p2_values)
    return final_probabilities

p1 = first['SurvivalProbability']

avg_x_data = []
avg_y_data = []
p2_values = []
final_probabilities = []

colors = plt.cm.viridis(np.linspace(0, 1, len(migration_rates)))  # Generate colors

# Load data and compute averages
avg_variances = load_and_compute_averages(len(migration_rates))
print(avg_variances)
# Compute the probabilities
final_probabilities = compute_probabilities(avg_variances, p1, slope, intercept)

# for idx, migration_rate in enumerate(migration_rates):
#     # Load data
#     peak_time_path = f"peak_time_difference_migration_rate_idx_{idx + 2}.csv"
#     peak_time_data = pd.read_csv(peak_time_path)['PeakTimeDifference'].values

#     variance_diff_path = f"variance_difference_migration_rate_idx_{idx + 2}.csv"
#     variance_diff_data = pd.read_csv(variance_diff_path)['VarianceDifference'].values

#     avg_y_data.append(np.mean(variance_diff_data))

#     # Calculate p2
#     p2 = p1 + slope * np.mean(variance_diff_data)
#     p2_values.append(p2)

#     # Calculate final probability
#     final_prob = 1 - (1 - p1) * (1 - p2)
#     final_probabilities.append(final_prob)


print(p2_values, p1)
ax.plot(migration_rates, final_probabilities, 'o-', label='test' )

ax.set_xscale('log')
ax.legend(frameon=False, loc=(0.1, 0.01))
ax.set_xlim(1e-8, 1)
ax.set_ylim(bottom=0)
ax.set_xlabel(r'migration rate (units: $\gamma$)')
ax.set_ylabel('survival probability')

#######################
plt.tight_layout()
plt.savefig('twodeme_trajectories_and_probs_with_theory.pdf', format='pdf', dpi=300)
plt.show()
print('Done')