using Serialization
using Plots

MIGRATION_RATES = exp10.(LinRange(-6,-0.5,10))

# Load the serialized results
results_file = "plotting_data/analysis_results.jld2"
results = open(deserialize, results_file)

# Prepare data for plotting
plot_data = Dict()
for (key, value) in results
    adjacency_matrix_idx, migration_rate_idx = key
    survival_probability, total_survivals, num_replicates = value
    standard_error = sqrt(survival_probability * (1 - survival_probability) / num_replicates)
    
    if !haskey(plot_data, adjacency_matrix_idx)
        plot_data[adjacency_matrix_idx] = ([], [], [])
    end
    push!(plot_data[adjacency_matrix_idx][1], MIGRATION_RATES[migration_rate_idx])
    push!(plot_data[adjacency_matrix_idx][2], survival_probability)
    push!(plot_data[adjacency_matrix_idx][3], standard_error)
end

# Plotting
group = [1, 3, 2, 2, 2, 1, 2, 2, 1, 1, 2]
group_colors = [:red, :green, :purple]
p = plot(legend=:topright, xscale=:log10)
for (adjacency_matrix_idx, (rates, probabilities, errors)) in plot_data
    sort_idx = sortperm(rates)  # To ensure the rates are in ascending order
    plot!(p, rates[sort_idx], probabilities[sort_idx], ribbon=errors[sort_idx], label="Network $(adjacency_matrix_idx)", color=group_colors[group[adjacency_matrix_idx]])
end

xlabel!(p, "Migration Rate")
ylabel!(p, "Probability of Survival")
title!(p, "Survival Probability vs. Migration Rate")

# Display or save the plot
display(p)
savefig(p, "plotting_data/survival_probability_plot.png")  # Uncomment to save the plot
