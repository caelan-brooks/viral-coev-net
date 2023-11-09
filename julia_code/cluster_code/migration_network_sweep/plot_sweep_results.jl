using Plots
using Serialization

# Load the preprocessed data
results_path = "./migration_network_sweep_results/migration_network_sweep_results.jld2"
data = open(deserialize, results_path)

survival_probabilities = data["survival_probabilities"]
num_replicates = data["num_replicates"]
network_sizes = 2:6
migration_rates = exp10.(LinRange(-6, 0.5, 10))

# Function to calculate standard error
function standard_error(p, n)
    se = sqrt.(p .* (1 .- p) ./ n)
    return se
end

# Initialize the plot with a transparent legend
p = plot(xscale=:log10, yscale=:lin, xlabel="Migration Rate", ylabel="Survival Probability", legend=:bottomleft, legendalpha=0.5)

# Plot the survival probability vs migration rate for each network size with error bars
for (idx, network_size) in enumerate(network_sizes)
    p_values = survival_probabilities[idx, :]
    errors = standard_error(p_values, num_replicates[idx, :])
    plot!(p, migration_rates, p_values, ribbon=errors, label="Network Size $network_size", marker=:circle)
end

# Add horizontal dashed line at 0.0904
hline!(p, [0.0904], linestyle=:dash, color=:black, label="Threshold")

# Display the plot
display(p)

# Save the figure
savefig(p, "survival_probability_vs_migration_rate.png")
