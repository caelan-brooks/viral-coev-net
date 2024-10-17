import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import matplotlib.patches as mpatches
from scipy.integrate import quad

# Set the global font to Times New Roman
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 20  # Default font size
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 18

def ridgeline(data, overlap=0.0, migration_rates=None, n_bins=500, yscale=50, hist_scale=0.002):
    if overlap > 1 or overlap < 0:
        raise ValueError('overlap must be in [0, 1]')
    ys = []

    # Filter out entries with zero antigenic variance from both deme1 and deme2
    filtered_data = [(list(filter(lambda x: x != 0, deme1)), list(filter(lambda x: x != 0, deme2)))
                     for deme1, deme2 in data]

    # Determine the range for bins
    data_min = min([min(deme) for deme1, deme2 in filtered_data for deme in [deme1, deme2] if deme])
    data_max = max([max(deme) for deme1, deme2 in filtered_data for deme in [deme1, deme2] if deme])
    bins = np.linspace(data_min, data_max, n_bins)
    bins = np.sort(np.unique(np.append(bins, 0)))  # Ensure zero is a bin edge
    
    for i, (deme1, deme2) in enumerate(filtered_data):
        y = i * (1.0 - overlap) * yscale
        ys.append(y)

        # Process each deme for histogram
        for deme, color in zip((deme1, deme2), ('darkgreen', 'purple')):
            hist, _ = np.histogram(deme, bins=bins, density=False)
            hist = hist.astype(float) * hist_scale
            plt.fill_between(bins[:-1], hist + y, y, step='post', alpha=0.3, color=color)
            plt.step(bins[:-1], hist + y, c=color, where='post')

        plt.xlim(0.1, 0.35)  # Adjusted xlim to focus on non-zero variance

    # Set y-axis labels if migration rates are provided
    if migration_rates is not None:
        plt.yticks(ys, ["{:.1e}".format(rate) for rate in migration_rates])

    plt.xlabel('Antigenic Diversity')  # Add x-axis label
    plt.ylabel('Migration Rate')  # Add y-axis label
    plt.title('Ridgeline Plot of Antigenic Diversity by Migration Rate')  # Add plot title


# Parameters
num_migration_rates = 11  # Adjust based on your data
csv_output_directory = "path/to/your/csv"  # Update with actual path

# Load data and prepare for plotting
data = []
labels = []
for migration_rate_idx in range(2, num_migration_rates+2):
    deme1_data = pd.read_csv(f"antigenic_variance_deme1_migration_rate_idx_{migration_rate_idx}.csv")['AntigenicVariance'].values
    deme2_data = pd.read_csv(f"antigenic_variance_deme2_migration_rate_idx_{migration_rate_idx}.csv")['AntigenicVariance'].values
    data.append((deme1_data, deme2_data))
    labels.append(f'Migration Rate {migration_rate_idx}')

# Create ridgeline plot
plt.figure(figsize=(12, 8))
migration_rates = 10**(np.linspace(-7, -0.5, num_migration_rates))
ridgeline(data, overlap=0.0, migration_rates=migration_rates, yscale=10, hist_scale=0.01)

plt.xlabel('Antigenic Diversity at Outbreak Peak')
plt.ylabel('Migration Rate')  # If you have a ylabel
plt.title('Distribution of Antigenic Diversity Across Migration Rates and Demes')

# Create custom legend
green_patch = mpatches.Patch(color='darkgreen', label='Deme 1', alpha=0.3)
purple_patch = mpatches.Patch(color='purple', label='Deme 2', alpha=0.3)
plt.legend(handles=[green_patch, purple_patch],loc='upper right')

plt.savefig('ridgeline_variance.pdf', format='pdf', dpi=300)
plt.savefig('ridgeline_variance.png', format='png')

# Constants
D = 0.01
N0 = 100
Nh = 2 * 10**6
R0 = 1.15
r = 1.5

# Function to calculate the integral for T
def integrand(s, k, N0, r, R0):
    return R0**(-(-1 + np.exp(r * s)) * k * N0 / r)

# Function to calculate T for a given migration rate k
def calculate_T(k, N0, r, R0):
    result, _ = quad(integrand, 0, 90, args=(k, N0, r, R0))
    return result

# Assuming 'data' and 'migration_rates' are already loaded as in your script

# Calculate average antigenic variance for each deme and migration rate
average_variances_deme1 = [np.mean(deme1) for deme1, _ in data]
average_variances_deme2 = [np.mean(deme2) for _, deme2 in data]

# Calculate the differences (deme2 - deme1) for each migration rate
differences = np.array(average_variances_deme2) - np.array(average_variances_deme1)

# Calculate new theoretical estimates
theoretical_estimates = np.array([2 * D * calculate_T(k, N0, r, R0) for k in migration_rates])

# Plotting the differences and new theoretical estimates
plt.figure(figsize=(12, 8))
plt.loglog(migration_rates, differences, marker='o', linestyle='-', color='blue', label='Empirical Difference')
plt.loglog(migration_rates, theoretical_estimates, linestyle='--', color='red', label='New Theoretical Estimate')

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Difference in Average Antigenic Variance (Deme 2 - Deme 1)')
plt.title('Difference in Average Antigenic Variance vs. Migration Rate')
plt.legend()
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('average_variance_difference_with_new_theory_loglog.pdf', format='pdf', dpi=300)
plt.savefig('average_variance_difference_with_new_theory_loglog.png', format='png')


# Calculate T for each migration rate again (for clarity and completeness)
Ts = np.array([calculate_T(k, N0, r, R0) for k in migration_rates])

# Calculate the ratio using N0 exp(rT)
ratios = N0 * np.exp(r * Ts)

# Plotting the ratio as a function of the migration rate on a log-log scale
plt.figure(figsize=(12, 8))
plt.loglog(migration_rates,migration_rates * ratios, marker='o', linestyle='-', color='green', label='Theoretical Ratio (Deme 1 / Deme 2)')

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Theoretical Rate of Diversity Flux')
plt.title('Theoretical Rate of Diversity Flux vs. Migration Rate')
plt.legend()
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('theoretical_ratio_infected_numbers_loglog.pdf', format='pdf', dpi=300)
plt.savefig('theoretical_ratio_infected_numbers_loglog.png', format='png')

# Calculate standard deviation (or standard error) for each deme and migration rate
std_devs_deme1 = [np.std(deme1, ddof=0) for deme1, _ in data]
std_devs_deme2 = [np.std(deme2, ddof=0) for _, deme2 in data]

# Calculate the standard error for the difference as sqrt(se1^2 + se2^2)
# Assuming independence between deme1 and deme2 variances
errors = np.sqrt(np.array(std_devs_deme1)**2 + np.array(std_devs_deme2)**2)

# Plotting the differences with error bars on a log-log scale
plt.figure(figsize=(12, 8))
plt.errorbar(migration_rates, differences, yerr=errors, fmt='o', linestyle='-', color='blue', label='Empirical Difference', capsize=5)

# Include the theoretical estimates if needed
plt.plot(migration_rates, theoretical_estimates, linestyle='--', color='red', label='New Theoretical Estimate')

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Difference in Average Antigenic Variance (Deme 2 - Deme 1)')
plt.title('Difference in Average Antigenic Variance vs. Migration Rate')
plt.legend()
plt.xscale("log")
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('variance_difference_loglog_with_errors.pdf', format='pdf', dpi=300)
plt.savefig('variance_difference_loglog_with_errors.png', format='png')


# Calculate the number of measurements for each migration rate
num_measurements_deme1 = [len(deme1) for deme1, _ in data]
num_measurements_deme2 = [len(deme2) for _, deme2 in data]
total_measurements = np.array(num_measurements_deme2)

# Plotting the total number of measurements for each migration rate
plt.figure(figsize=(12, 8))
plt.plot(migration_rates, total_measurements/num_measurements_deme1, marker='o', linestyle='-', color='purple', label='Total Measurements')
kstar = 20 * 1 / (2.5 * 1.37 * 10**6 * 2)
plt.plot(migration_rates, 1 - np.exp(-migration_rates/kstar))

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Number of Variance Measurements')
plt.xscale("log")
plt.yscale("linear")
plt.title('Number of Variance Measurements vs. Migration Rate')
plt.legend()
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('number_of_variance_measurements.pdf', format='pdf', dpi=300)
plt.savefig('number_of_variance_measurements.png', format='png')


plt.figure(figsize=(12, 8))
plt.plot(migration_rates, Ts, marker='o', linestyle='-', color='darkorange', label='Time T')

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Time T')
plt.title('Time T vs. Migration Rate')
plt.xscale('log')  # Setting x-axis to log scale if migration rates span several orders of magnitude
plt.yscale('log')  # Setting y-axis to log scale if T values span several orders of magnitude
plt.legend()
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('time_T_vs_migration_rate_loglog.pdf', format='pdf', dpi=300)
plt.savefig('time_T_vs_migration_rate_loglog.png', format='png')

# Theoretical prediction for deme 1
theoretical_deme1 = 2 * D * (np.log(Nh/N0) / r)
# Assuming Ts contains the calculated T values for each migration rate
theoretical_deme2 = theoretical_deme1 + 2 * D * Ts

plt.figure(figsize=(12, 8))

# Plotting average antigenic variance for deme 1 and deme 2
plt.plot(migration_rates, average_variances_deme1, marker='o', linestyle='-', color='blue', label='Average Variance Deme 1')
plt.plot(migration_rates, average_variances_deme2, marker='o', linestyle='-', color='green', label='Average Variance Deme 2')

# Plotting theoretical predictions
plt.axhline(y=theoretical_deme1, color='red', linestyle='--', label='Theoretical Deme 1')
plt.plot(migration_rates, theoretical_deme2, color='purple', linestyle='--', label='Theoretical Deme 2')

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Antigenic Variance')
plt.title('Antigenic Variance vs. Migration Rate')
plt.xscale('log')  # Applying log scale to x-axis if necessary
plt.yscale('linear')  # Applying log scale to y-axis if the range of values is wide
plt.legend()
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('average_variance_with_theory.pdf', format='pdf', dpi=300)
plt.savefig('average_variance_with_theory.png', format='png')

# Assuming the original migration_rates array is defined
min_rate = np.min(migration_rates)
max_rate = np.max(migration_rates)
fine_migration_rates = np.logspace(np.log10(min_rate), np.log10(max_rate), num=500)  # 500 points for finer resolution

# Recalculate T_analytical and T_asymptotic with the finer migration rate vector
T_analytical_fine = np.log(r / (fine_migration_rates * N0 * np.log(R0))) / r
T_asymptotic_fine = 1 / (fine_migration_rates * N0 * np.log(R0))

# Update theoretical predictions with the finer scale migration rate vector
variance_difference_theoretical_fine = 2 * D * T_analytical_fine
variance_difference_asymptotic_fine = 2 * D * T_asymptotic_fine

plt.figure(figsize=(12, 8))

# Plotting empirical variance difference with original migration rates
plt.errorbar(migration_rates, differences, yerr=errors, fmt='o', linestyle='-', color='blue', label='Empirical Difference', capsize=5)

# Plotting theoretical predictions with the finer scale migration rate vector
plt.plot(fine_migration_rates, variance_difference_theoretical_fine, linestyle='--', color='purple', label='New Theoretical Estimate')
plt.plot(fine_migration_rates, variance_difference_asymptotic_fine, linestyle='--', label='Asymptotic Theoretical Estimate')

# Keep the original theoretical estimates plotted against the original migration rates for comparison
plt.plot(migration_rates, theoretical_estimates, linestyle='--', color='red', label='Original Theoretical Estimate')

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Variance Difference (Deme 2 - Deme 1)')
plt.title('Variance Difference vs. Migration Rate with Theoretical Predictions')
plt.xscale('log')
plt.yscale('linear')
plt.ylim(0, 1.2 * max(differences))
plt.legend()
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('variance_difference_with_new_theory_fine.pdf', format='pdf', dpi=300)
plt.savefig('variance_difference_with_new_theory_fine.png', format='png')

# Calculate the variance of the antigenic variance for each deme at each migration rate
variances_deme1 = [np.var(deme1, ddof=1) for deme1, _ in data]  # Using ddof=1 for an unbiased estimator
variances_deme2 = [np.var(deme2, ddof=1) for _, deme2 in data]

# Assuming the calculation of variances_deme1 and variances_deme2 for the original migration rates
# Calculate the theoretical prediction for the finer migration rates
theoretical_variance_difference_fine = 4 * D**2 / (fine_migration_rates**2 * N0**2 * np.log(R0)**2)
theoretical_variance_difference = 4 * D**2 / (migration_rates**2 * N0**2 * np.log(R0)**2)

plt.figure(figsize=(12, 8))

# Plotting variance of antigenic variance for each deme with original migration rates
# plt.plot(migration_rates, variances_deme1, marker='o', linestyle='-', color='blue', label='Variance in Deme 1')
plt.plot(migration_rates, np.array(variances_deme2) - np.array(variances_deme1), marker='o', linestyle='-', color='green', label='Variance in Deme 2')
# plt.plot(migration_rates, np.array(variances_deme1) + theoretical_variance_difference, marker='o', linestyle='-', color='green', label='Variance in Deme 2')


# Plotting theoretical prediction for the difference in variance on the finer scale
plt.plot(fine_migration_rates, theoretical_variance_difference_fine, linestyle='--', color='red', label='Theoretical Variance Difference')

# Adding plot labels and title
plt.xlabel('Migration Rate')
plt.ylabel('Variance of Antigenic Variance')
plt.title('Variance of Antigenic Variance Across Samples vs. Migration Rate')
plt.xscale('log')  # Applying log scale to x-axis
plt.yscale('log')  # Applying log scale to y-axis due to the range of theoretical values
# plt.ylim(0,np.max(variances_deme2))
plt.legend()
plt.grid(True, which="both", ls="--")

# Save and show plot
plt.savefig('variance_of_antigenic_variance_with_theory.pdf', format='pdf', dpi=300)
plt.savefig('variance_of_antigenic_variance_with_theory.png', format='png')

def plot_average_peak_time_difference(migration_rates):
    num_migration_rates = len(migration_rates)
    average_peak_time_differences = []

    for migration_rate_idx in range(2, num_migration_rates + 2):
        # Load peak time differences data
        file_path = f"peak_time_difference_migration_rate_idx_{migration_rate_idx}.csv"
        peak_time_data = pd.read_csv(file_path)['PeakTimeDifference'].values
        
        # Calculate the average peak time difference
        average_peak_time = np.mean(peak_time_data)
        average_peak_time_differences.append(average_peak_time)
    
    # Plotting
    plt.figure(figsize=(12, 8))
    plt.plot(migration_rates, average_peak_time_differences, marker='o', linestyle='-', color='blue')
    plt.plot(migration_rates, 1/r * np.log(20 / 2 / N0 / migration_rates),linestyle='--')
    plt.xlabel('Migration Rate')
    plt.ylabel(r'Average Peak Time Difference, $\Delta$')
    plt.title('Average Peak Time Difference vs. Migration Rate')
    plt.xscale('log')  # Assuming migration rates span several orders of magnitude
    # plt.grid(True, which="both", linestyle='--', linewidth=0.5)
    
    plt.savefig('average_peak_time_difference.pdf', format='pdf', dpi=300)
    plt.savefig('average_peak_time_difference.png', format='png')
    
# Call the function
plot_average_peak_time_difference(migration_rates)

def plot_scatter_peak_time_vs_variance_diff(migration_rates):
    # Prepare a colormap and figure
    plt.figure(figsize=(12, 8))
    cmap = plt.get_cmap('viridis')
    colors = cmap(np.linspace(0, 1, len(migration_rates)))

    for idx, migration_rate in enumerate(migration_rates):
        # Load peak time differences data
        peak_time_path = f"peak_time_difference_migration_rate_idx_{idx + 2}.csv"
        peak_time_data = pd.read_csv(peak_time_path)['PeakTimeDifference'].values

        # Load variance differences data
        variance_diff_path = f"variance_difference_migration_rate_idx_{idx + 2}.csv"
        variance_diff_data = pd.read_csv(variance_diff_path)['VarianceDifference'].values

        # Calculate x data as 2*D*peak_time_difference
        x_data = 2 * D * peak_time_data

        # Plot scatter for each migration rate
        plt.scatter(x_data, variance_diff_data, color=colors[idx], alpha=0.5, label=f'Migration Rate {migration_rate}')

        # Calculate and plot averages
        avg_x = np.mean(x_data)
        avg_y = np.mean(variance_diff_data)
        plt.scatter(avg_x, avg_y, color=colors[idx], edgecolor='black', s=100, marker='o', label=f'Avg for {migration_rate}')

    # Plot y = x line
    max_limit = max(plt.xlim()[1], plt.ylim()[1])
    plt.plot([0, max_limit], [0, max_limit], 'k--', label='y = x')  # black dashed line

    plt.xlabel('2D deltaT')
    plt.ylabel('Variance Difference')
    plt.xlim(0,0.5)
    # plt.legend()
    plt.grid(True, which="both", linestyle='--', linewidth=0.5)

    plt.savefig('scatter_peak_time_vs_variance_diff.pdf', format='pdf', dpi=300)
    plt.savefig('scatter_peak_time_vs_variance_diff.png', format='png')

plot_scatter_peak_time_vs_variance_diff(migration_rates)

def plot_averages_only(migration_rates):
    # Prepare a colormap
    plt.figure(figsize=(12, 8))
    cmap = plt.get_cmap('viridis')
    colors = cmap(np.linspace(0, 1, len(migration_rates)))

    # Initialize lists for averages
    avg_x_data = []
    avg_y_data = []

    for idx, migration_rate in enumerate(migration_rates):
        # Load peak time differences data
        peak_time_path = f"peak_time_difference_migration_rate_idx_{idx + 2}.csv"
        peak_time_data = pd.read_csv(peak_time_path)['PeakTimeDifference'].values

        # Load variance differences data
        variance_diff_path = f"variance_difference_migration_rate_idx_{idx + 2}.csv"
        variance_diff_data = pd.read_csv(variance_diff_path)['VarianceDifference'].values

        # Calculate averages and append to lists
        avg_x_data.append(np.mean(2 * D * peak_time_data))
        avg_y_data.append(np.mean(variance_diff_data))

        # Label migration rates in scientific notation
        migration_rate_label = f"{migration_rate:.1e}"

        # Plot averages with simplified labels
        plt.scatter(avg_x_data[-1], avg_y_data[-1], color=colors[idx], edgecolor='black', s=100, label=migration_rate_label)

    xs = np.linspace(0, max(avg_x_data),200) / 2 / D
    ys = 2 * D * xs * (1 - np.exp(-20/2/N0 * (10 - xs)))
    plt.plot(2*D*xs, ys)


    # Plot y = x line
    max_limit = max(max(avg_x_data), max(avg_y_data))
    plt.plot([0, max_limit], [0, max_limit], 'k--', label='y = x')  # black dashed line

    plt.xlabel('2D * Average Peak Time Difference')
    plt.ylabel('Average Variance Difference')
    plt.title('Averages of 2D * Peak Time Difference vs. Variance Difference')
    plt.legend(title="Migration Rate")
    plt.grid(True, which="both", linestyle='--', linewidth=0.5)

    plt.savefig('averages_only.pdf', format='pdf', dpi=300)
    plt.savefig('averages_only.png', format='png')
# Call the function with the migration_rates array
plot_averages_only(migration_rates)
plt.show()
