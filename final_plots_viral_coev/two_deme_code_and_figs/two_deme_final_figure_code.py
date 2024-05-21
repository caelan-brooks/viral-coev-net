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

# First figure will have the trajectories and the survival probabilities
fig, axs = plt.subplots(1, 3, figsize=(8.5, 3.0))  # Adjust the figure size to balance the subplot shapes
labels = ['(a)', '(b)', '(c)']

for ax, label in zip(axs.flat, labels):
    ax.text(-0.01, 1.05, label, transform=ax.transAxes, va='bottom')

# Adjusting code for panel (a)
ax = axs[0]
trajectories_df = pd.read_csv("../trajectories_migration_rate_idx_1.csv")

# Initialize variables to store the maximum infected number and the corresponding time for the final trajectory
max_infected_number = 0
time_of_max_infection = 0

# Plot all trajectories
for index, row in trajectories_df.iterrows():
    trajectory = ast.literal_eval(row['parent'])  # Convert string representation of list to actual list
    time_points = np.linspace(0, 100, len(trajectory))
    ax.plot(time_points, trajectory, linewidth=1.5)
    
    # Update the maximum infected number and the time for the final trajectory
    if index == len(trajectories_df) - 1:
        max_infected_number = max(trajectory)
        time_of_max_infection = time_points[np.argmax(trajectory)]

print(time_of_max_infection)
# Add vertical line at the time of maximum infection for the final trajectory
ax.axvline(x=time_of_max_infection, color='black', linestyle='--', linewidth=1)
ax.text(time_of_max_infection, 0.2, '$T$', fontsize=written_text_fontsize, ha='center', va='bottom')

# Add text and arrow labeling the peak as the "outbreak peak"
ax.annotate('outbreak peak', xy=(time_of_max_infection, max_infected_number), xycoords='data',
            xytext=(time_of_max_infection + 7, max_infected_number), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.05, width=0.5, headwidth=3, headlength=3, connectionstyle='arc3,rad=0'),
            fontsize=written_text_fontsize, ha='left', va='center')


start_time = 65
end_time = 100

ax.text(x=start_time + (end_time - start_time) / 3.0, y=1e6, s="escape", rotation='horizontal', ha='center', va='center', color='green', fontsize=written_text_fontsize)
ax.text(x=end_time * 0.83, y=3, s="extinction", ha='center', va='center', color='red', fontsize=written_text_fontsize)

ax.set_yscale('log')
ax.set_xlabel(r'time (units: $\gamma^{-1}$)')
ax.set_ylabel(r'total infected number, $N$')
ax.set_xlim(0, 100)
ax.set_ylim(bottom=1)


# Adjusting code for panel (b)
ax = axs[2]  # This selects the last panel position in a 3x2 grid

# Load data from CSV
data = pd.read_csv("../variances_and_probabilities.csv")

# Extract variance and probability
variance = data['variance'].values
probability = data['probability'].values
error = data['error'].values

# Perform linear regression
A = np.vstack([variance, np.ones(len(variance))]).T
slope, intercept = np.linalg.lstsq(A, probability, rcond=None)[0]

# Define linear fit function
linear_fit = lambda x: slope * x + intercept

# Scatter plot with error bars
ax.errorbar(variance, probability, yerr=error, fmt='o', label="Survival Probability", color='blue', ecolor='lightgray', elinewidth=1, capsize=1)
print(error)

# Plot linear fit
x_vals = np.linspace(min(variance), max(variance), 100)
ax.plot(x_vals, linear_fit(x_vals), 'r--', label=f"Linear Fit\ny = {slope:.3f}x + {intercept:.3f}")

# Customize plot
ax.set_xlabel(r"antigenic diversity at outbreak peak, $V(T)$")
ax.set_ylabel("escape probability", labelpad=1)
ax.set_ylim(0, 0.5)

# Function to create a Gaussian plot
def plot_gaussian(ax, mean, variance, color):
    x = np.linspace(mean - 9 * np.sqrt(variance), mean + 9 * np.sqrt(variance), 100)
    y = (1 / (np.sqrt(2 * np.pi * variance))) * np.exp(-0.5 * ((x - mean) ** 2 / variance))
    ax.plot(x, y, color=color)
    ax.set_xlim(mean - 3 * 1.0, mean + 3 * 1.0)
    ax.set_ylim(0, max(y) * 1.1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("$x$")
    ax.set_ylabel("$n(x,T)$")
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    return max(y)

# Adding inset for low variance Gaussian
low_variance = variance[0] / 3
mean = 0
inset_ax_low = ax.inset_axes([0.15, 0.6, 0.3, 0.3])  # Adjust position and size as needed
low_ymax = plot_gaussian(inset_ax_low, mean, low_variance, 'black')

# Adding inset for high variance Gaussian
high_variance = variance[-1] * 2
inset_ax_high = ax.inset_axes([0.67, 0.15, 0.3, 0.3])  # Adjust position and size as needed
high_ymax = plot_gaussian(inset_ax_high, mean, high_variance, 'black')

# Set the same ylim for both insets
common_ylim = max(low_ymax, high_ymax) * 1.1
inset_ax_low.set_ylim(0, common_ylim)
inset_ax_high.set_ylim(0, common_ylim)

# Adding arrows pointing to data points
high_var_idx = np.argmax(variance)
low_var_idx = np.argmin(variance)
ax.annotate('', xy=(variance[high_var_idx], probability[high_var_idx]), xytext=(0.85, 0.25), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.1, width=0.5, headwidth=5, headlength=5))
ax.annotate('', xy=(variance[low_var_idx], probability[low_var_idx]), xytext=(0.25, 0.75), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.1, width=0.5, headwidth=5, headlength=5))

# Customize main plot
ax.set_xticks(np.linspace(min(variance), max(variance), num=6))  # Generate xticks
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:.2f}'))  # Format xticks with 1 decimal place

ax = axs[1]

# Parameters for the Gaussians
mean = 0
sigma1 = np.sqrt(1)  # Standard deviation for deme 1

# Function to calculate FWHM
def fwhm(sigma):
    return 2 * np.sqrt(2 * np.log(2)) * sigma

# Generate x values
x = np.linspace(-5, 5, 1000)

# Gaussian distributions
y1 = norm.pdf(x, mean, sigma1)

# Scale for visibility if needed, or adjust plot limits or aspect
scaling_factor = 7
y1_scaled = y1 * scaling_factor

# Plotting on the specified axis
ax.plot(x, y1_scaled, linewidth=2, color='blue')

# Shade under the curve
ax.fill_between(x, y1_scaled, color='blue', alpha=0.1)

# Calculate FWHM for each distribution
fwhm1 = fwhm(sigma1)

# Height for arrows (half-maximum of the unscaled distributions)
height1 = max(y1) * scaling_factor / 2

# Plot two-way arrows for FWHM
ax.annotate('', xy=(mean - fwhm1/2, height1), xytext=(mean + fwhm1/2, height1),
            arrowprops=dict(arrowstyle='<->', color='black'))


# Text annotations for FWHM
ax.text(mean, height1 * 0.95, r'$\sqrt{V(t)}$', ha='center', va='top', color='black')

# Annotation for the total number of infected hosts
# ax.annotate(r'area = $N(t)$', xy=(0, 0.5), xycoords='data',
#             xytext=(-3.5, height1 * 0.8), textcoords='data',
#             arrowprops=dict(facecolor='black', shrink=0.01, width=0.5, headwidth=3, headlength=3),
#             fontsize=written_text_fontsize, ha='center', va='center', color='black')
ax.text(mean, height1 * 0.25, r'area = $N(t)$', ha='center', va='top', color='black')

# Additional plot settings
ax.set_xlabel(r'antigenic coordinate, $x$')
ax.set_ylabel(r'infected number density $n(x,t)$')
ax.set_yticks([])
ax.set_xticks([])
ax.set_ylim(bottom=0)

plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.1) 
plt.subplots_adjust(right=0.96)
plt.savefig('single_deme.pdf', format='pdf', dpi=300)
###############################################################

migration_rates = 10**(np.linspace(-10, 1.0, 12))
N0 = 100
D = 1/100

# Read the CSV file
df = pd.read_csv("../analysis_results.csv")

# Define colors for horizontal lines and ribbons
color_first = 'royalblue'
color_last = 'darkorange'

# Define extended migration rate range for shading
migration_rate_min = 1e-12
migration_rate_max = 30

# Calculate the standard error for a binomial distribution
n = 10000  # number of replicates
df['StandardError'] = np.sqrt(df['SurvivalProbability'] * (1 - df['SurvivalProbability']) / n)

# Extracting the first, last, and middle elements
first = df.iloc[0]
last = df.iloc[-1]
middle_df = df.iloc[1:-1]

# New data
df_new = pd.read_csv("../analysis_results_noback.csv")
n_new = 10000  # Assuming the same number of replicates for the new data
df_new['StandardError'] = np.sqrt(df_new['SurvivalProbability'] * (1 - df_new['SurvivalProbability']) / n_new)
middle_df_new = df_new.iloc[1:-3]  # Adjust as per your new data structure


fig, axs = plt.subplots(1, 3, figsize=(8.5, 3.0))  # Adjust the figure size to balance the subplot shapes
labels = ['(a)', '(b)', '(c)']

for ax, label in zip(axs.flat, labels):
    # Position: x=0, y=1 in axis coordinates (top left corner), transform=ax.transAxes
    # uses the axes coordinate system (0,0 is bottom left of the plot and 1,1 is top right)
    ax.text(-0.02, 1.05, label, transform=ax.transAxes, va='bottom')


ax = axs[0]

# Plotting code adjusted for the subplot
# ax.errorbar(middle_df['MigrationRate'], middle_df['SurvivalProbability'], yerr=middle_df['StandardError'], fmt='o', capsize=5, color='saddlebrown', label='Two-way Migration (1 \u2194 2)')
# ax.errorbar(middle_df_new['MigrationRate'], middle_df_new['SurvivalProbability'], yerr=middle_df_new['StandardError'], fmt='v', capsize=5, color='darkgreen', label='One-way Migration (1 \u2192 2)')
ax.scatter(middle_df['MigrationRate'], middle_df['SurvivalProbability'], color='saddlebrown', label='Two-way Migration (1 \u2194 2)', marker='o')
ax.scatter(middle_df_new['MigrationRate'], middle_df_new['SurvivalProbability'], color='darkgreen', label='One-way Migration (1 \u2192 2)', marker='v')

ax.axhline(y=first['SurvivalProbability'], color='royalblue', linestyle='--')

ax.fill_between([migration_rate_min, migration_rate_max], first['SurvivalProbability'] - first['StandardError'], first['SurvivalProbability'] + first['StandardError'], color=color_first, alpha=0.3)

ax.axhline(y=last['SurvivalProbability'], color='darkorange', linestyle='--')
ax.fill_between([migration_rate_min, migration_rate_max], last['SurvivalProbability'] - last['StandardError'], last['SurvivalProbability'] + last['StandardError'], color=color_last, alpha=0.3)

def load_and_compute_averages(num_migration_rates):
    avg_variances = []
    for idx in range(2, num_migration_rates + 2):
        deme2_data = pd.read_csv(f"../antigenic_variance_deme2_migration_rate_idx_{idx}.csv")['AntigenicVariance'].values
        avg_variance = np.mean(deme2_data) * len(deme2_data) / 10000 + (-intercept/slope) * (1 - len(deme2_data)/10000) # adjust for the fact that at low migration not all trajectories result in secondary outbreaks
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
final_probabilities = compute_probabilities(avg_variances, p1, slope, intercept)

ax.plot(migration_rates, final_probabilities, 'k--', label='Theory' )

ax.set_xscale('log')
ax.legend(frameon=False, loc=(0.05, 0.7))
ax.set_xlim(migration_rate_min, migration_rate_max)
ax.set_ylim(bottom=0, top=0.7)
ax.set_xlabel(r'migration rate, $k$ (units: $\gamma$)')
ax.set_ylabel('escape probability')

ax = axs[1]

df = pd.read_csv("../single_trajectory.csv")

# Plot total infected for each deme
ax.plot(df['times'], df['total_infected_1'], color='blue', linewidth=2, label="deme 1")
ax.plot(df['times'], df['total_infected_2'], color='green', linewidth=2, label="deme 2")

max_idx_deme1 = df['total_infected_1'].idxmax()
max_idx_deme2 = df['total_infected_2'].idxmax()

time_max_infected_deme1 = df.loc[max_idx_deme1, 'times']
time_max_infected_deme2 = df.loc[max_idx_deme2, 'times']

# Vertical dashed lines
ax.axvline(time_max_infected_deme1, color='blue', linestyle='--', linewidth=2)
ax.axvline(time_max_infected_deme2, color='green', linestyle='--', linewidth=2)

# Draw a two-headed arrow between the lines
ax.annotate('', xy=(time_max_infected_deme1, ax.get_ylim()[1]*2), xytext=(time_max_infected_deme2, ax.get_ylim()[1]*2),
            arrowprops=dict(arrowstyle="<->", color='black'))

# Label the arrow as \Delta
mid_point = (time_max_infected_deme1 + time_max_infected_deme2) / 2
ax.text(mid_point, ax.get_ylim()[1]/0.3, r'$\Delta T$', horizontalalignment='center', color='black')
ax.text(time_max_infected_deme1, 0.2 , r'$T_1$', color='black', horizontalalignment='center')
ax.text(time_max_infected_deme2, 0.2 , r'$T_2$', color='black', horizontalalignment='center')

# Find the antigenic variances at the times of maximum infection
antigenic_variance_at_max1 = df.loc[max_idx_deme1, 'antigenic_variance_1']
antigenic_variance_at_max2 = df.loc[max_idx_deme2, 'antigenic_variance_2']

# Labeling
ax.set_ylabel(r'infected number, $N$')
ax.set_yscale('log')
ax.set_ylim(1, 10**(7.5))

# Add legends
ax.legend(loc='upper left', frameon=False)

# Labeling axes
ax.set_xlabel(r'time (units: $\gamma^{-1}$)')
ax.set_xlim(0,18)

# Create an inset figure for the Gaussians
inset_ax = ax.inset_axes([0.05, 0.03, 0.8, 0.45])  # Adjust the position and size of the inset as needed
inset_ax.patch.set_alpha(0.9)
# Generate x values for the Gaussians
x = np.linspace(-1.3, 1.3, 1000)

# Gaussian distributions with given variances
mean = 0
y1 = norm.pdf(x, mean, np.sqrt(antigenic_variance_at_max1))
y2 = norm.pdf(x, mean, 1 * np.sqrt(antigenic_variance_at_max2))

# Plotting the Gaussians on the inset axis
inset_ax.plot(x, y1, color='blue', linewidth=2, label=r'$n_1(x,T_1)$')
inset_ax.plot(x, y2, color='green', linewidth=2, label=r'$n_2(x,T_2)$')

# Remove ticks
inset_ax.set_xticks([])
inset_ax.set_yticks([])

# Set labels
# inset_ax.set_xlabel(r'$x$', fontsize=10)
# inset_ax.set_ylabel(r'$n_i(x,T_i)$', fontsize=10)
# Add labels inside the inset boundaries
inset_ax.text(0.5, 0.05, r'$x$', transform=inset_ax.transAxes, ha='center', va='center', fontsize=10)
inset_ax.text(0.05, 0.5, r'density, $n$', transform=inset_ax.transAxes, ha='center', va='center', rotation='vertical', fontsize=10)

# Add a legend
inset_ax.legend(loc='best', fontsize=8, frameon=False, handlelength=1)

ax = axs[2]
cmap = plt.get_cmap('inferno')
colors = cmap(np.linspace(0, 1, len(migration_rates)))

avg_x_data = []
avg_y_data = []

for idx, migration_rate in enumerate(migration_rates):
    # Adjust file paths as necessary
    peak_time_path = f"../peak_time_difference_migration_rate_idx_{idx + 2}.csv"
    peak_time_data = pd.read_csv(peak_time_path)['PeakTimeDifference'].values

    variance_diff_path = f"../variance_difference_migration_rate_idx_{idx + 2}.csv"
    variance_diff_data = pd.read_csv(variance_diff_path)['VarianceDifference'].values

    avg_x_data.append(np.mean(2 * D * peak_time_data))
    avg_y_data.append(np.mean(variance_diff_data))
    # print(np.mean(variance_diff_data))

    migration_rate_label = f"{migration_rate:.1e}"
    ax.scatter(avg_x_data[-1], avg_y_data[-1], color=colors[idx], edgecolor='black', s=100, label=migration_rate_label)

# Plot theoretical curve (adjust as needed based on your theory or model)
xs = np.linspace(0, max(avg_x_data), 200) / 2 / D
ys = 2 * D * xs * (1 - np.exp(-20/2/N0 * (10 - xs)))
# ax.plot(2 * D * xs, ys)

# Plot y = x line
max_limit = max(max(avg_x_data), max(avg_y_data))
ax.plot([0, max_limit], [0, max_limit], 'k--', label='y = x')  # Black dashed line

ax.set_xlabel(r'$\langle 2 D \Delta T \rangle$')
ax.set_ylabel(r'$\langle V_2(T_2)  - V_1(T_1) \rangle$')
# ax.set_title('Averages of 2D * Peak Time Difference vs. Variance Difference')
# ax.legend(title="Migration Rate")
ax.grid(True, which="both", linestyle='--', linewidth=0.5)

plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.05) 
# plt.tight_layout()
plt.subplots_adjust(right=0.99)
plt.savefig('two_deme.pdf', format='pdf', dpi=300)