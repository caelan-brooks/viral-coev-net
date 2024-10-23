import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
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

rcParams['figure.dpi'] = 100  # or set it to any constant value

# Create a directory for split figures if it doesn't exist
output_dir = 'split_figures'
os.makedirs(output_dir, exist_ok=True)

# Function to create and save the first subplot as its own figure
def plot_first_subplot():
    fig, ax = plt.subplots()
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

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'first_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_first_subplot()

# Function to create and save the second subplot as its own figure
def plot_second_subplot():
    fig, ax = plt.subplots()

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

    CROSS_REACTIVE_R = 3.0

    # Scatter plot with error bars
    ax.errorbar(variance / CROSS_REACTIVE_R**2, probability, yerr=error, fmt='o', label=None, color='blue', ecolor='lightgray', elinewidth=1, capsize=1)
    print(error)

    # Plot linear fit
    x_vals = np.linspace(min(variance/CROSS_REACTIVE_R**2), max(variance/CROSS_REACTIVE_R**2), 100)
    formatted_slope = f"{slope * CROSS_REACTIVE_R**2:.2f}"  # Format the slope with 2 decimal places
    ax.plot(x_vals, linear_fit(x_vals * CROSS_REACTIVE_R**2), 'r--', label=rf"slope = {formatted_slope}")

    # Customize plot
    ax.set_xlabel(r"diversity at outbreak peak, $V(T) / r_0^2$")
    ax.set_ylabel("escape probability", labelpad=1)
    ax.set_ylim(0, 0.5)
    ax.legend(frameon=False)

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
    ax.annotate('', xy=(variance[high_var_idx]/CROSS_REACTIVE_R**2, probability[high_var_idx]), xytext=(0.85, 0.25), textcoords='axes fraction',
                arrowprops=dict(facecolor='black', shrink=0.1, width=0.5, headwidth=5, headlength=5))
    ax.annotate('', xy=(variance[low_var_idx]/CROSS_REACTIVE_R**2, probability[low_var_idx]), xytext=(0.25, 0.75), textcoords='axes fraction',
                arrowprops=dict(facecolor='black', shrink=0.1, width=0.5, headwidth=5, headlength=5))

    # Customize main plot
    ax.set_xticks(np.linspace(min(variance/CROSS_REACTIVE_R**2), max(variance/CROSS_REACTIVE_R**2), num=4))  # Generate xticks
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:.3f}'))  # Format xticks with 1 decimal place

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'second_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_second_subplot()

# Function to create and save the third subplot as its own figure
def plot_third_subplot():
    fig, ax = plt.subplots()

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
    ax.text(mean, height1 * 0.25, r'area = $N(t)$', ha='center', va='top', color='black')

    # Additional plot settings
    ax.set_xlabel(r'antigenic coordinate, $x$')
    ax.set_ylabel(r'infected number density $n(x,t)$')
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylim(bottom=0)

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'third_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_third_subplot()

# Function to create and save the fourth subplot as its own figure
def plot_fourth_subplot():
    fig, ax = plt.subplots()

    # Load data from CSV
    data = pd.read_csv("../pop_size_scaling.csv")
    
    # Scatter plot of the data
    ax.scatter(data.HostPopulationSize, data.SurvivalProbability, label=None, color='blue')
    
    # Set x-axis to log scale
    ax.set_xscale("log")
    ax.set_xlabel(r"number of hosts, $N_h$")
    ax.set_ylabel("escape probability")
    ax.grid()

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'fourth_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_fourth_subplot()

def plot_second_figure_first_subplot():
    fig, ax = plt.subplots()

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
    df_new['StandardError'] = np.sqrt(df_new['SurvivalProbability'] * (1 - df_new['SurvivalProbability']) / n)
    middle_df_new = df_new.iloc[1:-3]

    # Scatter plots
    ax.scatter(middle_df['MigrationRate'], middle_df['SurvivalProbability'], color='saddlebrown', label='two-way migration (1 ↔ 2)', marker='o')
    ax.scatter(middle_df_new['MigrationRate'], middle_df_new['SurvivalProbability'], color='darkgreen', label='one-way migration (1 → 2)', marker='v')

    # Horizontal lines and shaded regions
    ax.axhline(y=first['SurvivalProbability'], color=color_first, linestyle='--')
    ax.fill_between([migration_rate_min, migration_rate_max], first['SurvivalProbability'] - first['StandardError'], first['SurvivalProbability'] + first['StandardError'], color=color_first, alpha=0.3)

    ax.axhline(y=last['SurvivalProbability'], color=color_last, linestyle='--')
    ax.fill_between([migration_rate_min, migration_rate_max], last['SurvivalProbability'] - last['StandardError'], last['SurvivalProbability'] + last['StandardError'], color=color_last, alpha=0.3)

    # Load data and compute averages
    def load_and_compute_averages(num_migration_rates):
        avg_variance_differences = []
        probability_of_spreading = []
        for idx in range(2, num_migration_rates + 2):
            deme2_data = pd.read_csv(f"../variance_difference_migration_rate_idx_{idx}.csv")['VarianceDifference'].values
            avg_variance_difference = np.mean(deme2_data)
            avg_variance_differences.append(avg_variance_difference)
            probability_of_spreading.append(len(deme2_data)/10000)
        return np.array(avg_variance_differences), np.array(probability_of_spreading)

    # Calculate slope using linear regression, similar to the second subplot code
    data_variance = pd.read_csv("../variances_and_probabilities.csv")
    variance = data_variance['variance'].values
    probability = data_variance['probability'].values
    A = np.vstack([variance, np.ones(len(variance))]).T
    slope, intercept = np.linalg.lstsq(A, probability, rcond=None)[0]

    # Calculate probabilities
    def compute_probabilities(avg_variance_differences, probability_of_spreading, p1, slope):
        p2_values = probability_of_spreading * (p1 + slope * avg_variance_differences)
        final_probabilities = 1 - (1 - p1) * (1 - p2_values)
        return final_probabilities

    p1 = first['SurvivalProbability']
    avg_variance_differences, probability_of_spreading = load_and_compute_averages(len(migration_rates))
    final_probabilities = compute_probabilities(avg_variance_differences, probability_of_spreading, p1, slope)

    # Plot the final theoretical probabilities
    ax.plot(migration_rates, final_probabilities, 'k--', label='theory')

    # Customize the plot
    ax.set_xscale('log')
    ax.legend(frameon=False, loc=(0.05, 0.7))
    ax.set_xlim(migration_rate_min, migration_rate_max)
    ax.set_ylim(bottom=0, top=0.7)
    ax.set_xlabel(r'migration rate, $k / \gamma$')
    ax.set_ylabel('escape probability')

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'second_figure_first_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_second_figure_first_subplot()

# Function to create and save the second subplot of the second figure as its own figure
def plot_second_figure_second_subplot():
    fig, ax = plt.subplots()

    # Load data from CSV
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

    print("t1 = ", time_max_infected_deme1)

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
    inset_ax.text(0.5, 0.05, r'$x$', transform=inset_ax.transAxes, ha='center', va='center', fontsize=10)
    inset_ax.text(0.05, 0.5, r'density, $n$', transform=inset_ax.transAxes, ha='center', va='center', rotation='vertical', fontsize=10)

    # Add a legend
    inset_ax.legend(loc='best', fontsize=8, frameon=False, handlelength=1)

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'second_figure_second_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_second_figure_second_subplot()

from matplotlib.colors import LogNorm

# Function to create and save the third subplot of the second figure as its own figure
def plot_second_figure_third_subplot():
    fig, ax = plt.subplots()
    r0 = 3

    # Load necessary data
    migration_rates = 10**(np.linspace(-10, 1.0, 12))
    D = 0.01
    cmap = plt.get_cmap('inferno')
    color_norm = LogNorm(vmin=min(migration_rates), vmax=max(migration_rates))
    colors = cmap(color_norm(migration_rates))  # Get colors based on log scale

    avg_x_data = []
    avg_y_data = []

    # Iterate over migration rates to load and plot data
    for idx, migration_rate in enumerate(migration_rates):
        peak_time_path = f"../peak_time_difference_migration_rate_idx_{idx + 2}.csv"
        peak_time_data = pd.read_csv(peak_time_path)['PeakTimeDifference'].values

        variance_diff_path = f"../variance_difference_migration_rate_idx_{idx + 2}.csv"
        variance_diff_data = pd.read_csv(variance_diff_path)['VarianceDifference'].values

        avg_x_data.append(np.mean(2 * D * peak_time_data))
        avg_y_data.append(np.mean(variance_diff_data))

        # Scatter plot points with migration rate-based colors
        ax.scatter(avg_x_data[-1] / r0**2, avg_y_data[-1] / r0**2, color=colors[idx], edgecolor='black', s=100, label=None)

    # Create the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=color_norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(r'$k/\gamma$', labelpad=0, rotation=0)
    cbar.ax.xaxis.set_label_position('top')

    # Set plot limits
    ax.set_ylim(bottom=0, top=0.25 / r0**2)

    # Define constants and functions for the theoretical model
    N0 = 100
    sig = 2
    beta = 2.5
    gamma = 1.0
    F = beta - gamma
    T1 = 10

    # Define the function cdf_delta
    def cdf_delta(t, k):
        return np.exp(-2 * N0 * k / sig**2 * (np.exp(F * t) - 1))

    # Time and ks values
    ts = np.linspace(0, 15, int(1e4))
    ks = np.logspace(-10, 0.5, 30)

    # Initialize avg_deltas array
    avg_deltas = np.full_like(ks, np.nan)

    # Calculate avg_deltas
    for ii in range(len(ks)):
        avg_delta = np.trapz(cdf_delta(ts, ks[ii]), ts)
        avg_deltas[ii] = avg_delta

    # Initialize var_diffs array
    var_diffs = np.full_like(ks, np.nan)

    # Calculate var_diffs
    for ii in range(len(ks)):
        term = -ks[ii] * (T1 - avg_deltas[ii]) * N0 * np.exp(F * avg_deltas[ii]) - sig**2 / 2 / beta
        var_diffs[ii] = 2 * D * avg_deltas[ii] * (1 - np.exp(term))

    # Plot y = x line
    max_limit = max(max(avg_x_data), max(avg_y_data))
    ax.plot([0, max_limit / r0**2], [0, max_limit / r0**2], 'k--', label='y = x')

    # Customize the plot
    ax.set_xlabel(r'$\langle 2 D \Delta T \rangle / r_0^2$')
    ax.set_ylabel(r'$\langle V_2(T_2)  - V_1(T_1) \rangle / r_0^2$')
    ax.grid(True, which="both", linestyle='--', linewidth=0.5)

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'second_figure_third_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_second_figure_third_subplot()

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

# Function to create and save the fourth subplot of the second figure as its own figure
def plot_second_figure_fourth_subplot():
    fig, ax = plt.subplots()

    # Function to calculate standard error
    def standard_error(p, n=10000):
        return np.sqrt(p * (1 - p) / n)

    # Load the data
    csv_file = "../survival_probabilities.csv"  # Adjust the file name/path if necessary
    data = pd.read_csv(csv_file)
    mutation_rates = np.linspace(0.001, 0.02, 10)
    migration_rates = np.power(10.0, np.linspace(-7, -0.5, 10))
    cmap = plt.get_cmap('inferno')
    color_norm = LogNorm(vmin=min(migration_rates), vmax=max(migration_rates))

    # Plot survival probability for different migration rates
    for migration_rate_idx in range(1, len(migration_rates) + 1):
        subset = data[data["MigrationRateIdx"] == migration_rate_idx]
        errors = subset["SurvivalProbability"].apply(lambda p: standard_error(p))
        ax.plot(mutation_rates, subset["SurvivalProbability"], '-o', color=cmap(color_norm(migration_rates))[migration_rate_idx-1], label=f'{migration_rates[migration_rate_idx - 1]:.1e}', markersize=4)

    ax.set_xlabel(r"mutation rate, $D$")
    ax.set_ylabel("escape probability")
    ax.set_xscale("linear")
    ax.grid(True)
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))

    # Create inset plot using ax.inset_axes
    ax_inset = ax.inset_axes([0.147, 0.49, 0.27, 0.5])  # Position and size of inset [x0, y0, width, height]
    ax_inset.patch.set_alpha(0.5)  # Set background transparent

    # Find the migration rate with the highest survival probability, skipping cases with any survival probability of 0 or 1
    best_migration_rates = []
    for idx, mutation_rate in enumerate(mutation_rates):
        subset = data[data["MutationRateIdx"] == idx + 1]  # Filter rows with the correct MutationRateIdx
        
        # Check if any survival probabilities are 0 or 1, if so, skip this mutation rate
        if (subset["SurvivalProbability"] <= 1e-3).any() or (subset["SurvivalProbability"] == 1).any():
            best_migration_rates.append(np.nan)  # Skip this mutation rate
            continue

        # Otherwise, find the max survival probability
        if not subset.empty:
            best_idx = subset["SurvivalProbability"].idxmax()  # Find index of max survival probability
            best_migration_idx = data.loc[best_idx, "MigrationRateIdx"]  # Get the MigrationRateIdx
            best_migration = migration_rates[best_migration_idx - 1]  # Look up the actual migration rate (adjust for 0-based index)
            best_migration_rates.append(best_migration)
        else:
            best_migration_rates.append(np.nan)  # Handle cases where no valid survival probabilities exist

    # Plot in the inset
    ax_inset.plot(mutation_rates, best_migration_rates, '-o', color='black', markersize=4)
    ax_inset.set_xlabel(r"$D$", fontsize=9)
    ax_inset.set_ylabel(r"optimal $k/\gamma$", fontsize=9, labelpad=0.0)  # Adjust the labelpad to reduce space
    ax_inset.set_xscale("linear")
    ax_inset.set_yscale("log")
    ax_inset.grid(True)

    # Set y-axis limits to reduce padding
    ax_inset.set_ylim(1.2e-6, 8e-5)
    ax_inset.set_xlim(min(mutation_rates[3:7]), max(mutation_rates[3:7]))  # Set x-axis limits for the inset to match mutation rates range

    # Adjust tick parameters to make ticks face inward
    def format_ticks(x, _):
        return f'{x:.2g}'

    # Set custom ticks and formatter for inset
    ax_inset.set_xticks([0.005, 0.01, 0.015])  # Custom tick locations
    ax_inset.xaxis.set_major_formatter(FuncFormatter(format_ticks))
    ax_inset.tick_params(axis='both', which='major', labelsize=6, direction='in')

    # Save the figure in the split_figures directory
    output_path_base = os.path.join(output_dir, 'second_figure_fourth_subplot')
    fig.savefig(f"{output_path_base}.png")
    fig.savefig(f"{output_path_base}.pdf")
    fig.savefig(f"{output_path_base}.svg")
    plt.close(fig)

# Example usage
plot_second_figure_fourth_subplot()
