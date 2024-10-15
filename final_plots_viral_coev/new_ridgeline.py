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

def ridgeline(data, overlap=0.0, migration_rates=None, n_bins=500, yscale=50, hist_scale=0.002, density=False):
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
            hist, _ = np.histogram(deme, bins=bins, density=density)
            hist = hist.astype(float) * hist_scale
            plt.fill_between(bins[:-1], hist + y, y, step='post', alpha=0.3, color=color)
            plt.step(bins[:-1], hist + y, c=color, where='post')

        plt.xlim(0.1, 0.35)  # Adjusted xlim to focus on non-zero variance

    # Set y-axis labels if migration rates are provided
    if migration_rates is not None:
        plt.yticks(ys, ["{:.1e}".format(rate) for rate in migration_rates])

    # plt.xlabel('Antigenic Diversity')  # Add x-axis label
    # plt.ylabel('Migration Rate')  # Add y-axis label
    # plt.title('Ridgeline Plot of Antigenic Diversity by Migration Rate')  # Add plot title


# Parameters
num_migration_rates = 12  # Adjust based on your data
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
migration_rates = 10**(np.linspace(-10, 1.0, num_migration_rates))
ridgeline(data, overlap=0.0, migration_rates=migration_rates, yscale=60, hist_scale=0.03)
# ridgeline(data, overlap=0.0, migration_rates=migration_rates, yscale=50, hist_scale=1.0)

plt.xlabel('antigenic diversity at outbreak peak')
plt.ylabel(r'migration rate, $k/\gamma$')  # If you have a ylabel
# plt.title('Distribution of Antigenic Diversity Across Migration Rates and Demes')

# Create custom legend
green_patch = mpatches.Patch(color='darkgreen', label='deme 1', alpha=0.3)
purple_patch = mpatches.Patch(color='purple', label='deme 2', alpha=0.3)
plt.legend(handles=[green_patch, purple_patch],loc='upper right')

plt.savefig('ridgeline_variance.pdf', format='pdf', dpi=300)
plt.savefig('ridgeline_variance.png', format='png')
