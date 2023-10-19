using Serialization
using Plots

const BETA_VALUES = LinRange(1.05, 6, 10)  # Adjust as needed
const DATA_PATH = "processed_data_betasweep/extinction_results.jld2"  # Adjust path as necessary

# Load preprocessed data
results = deserialize(open(DATA_PATH, "r"))

# Extract extinction probabilities and compute standard errors
n_files = [sum(values(results[idx])) for idx in 1:length(BETA_VALUES)]
demographic_extinctions = [results[idx]["demographic"] / n_files[idx] for idx in 1:length(BETA_VALUES)]
immune_extinctions = [results[idx]["immune"] / n_files[idx] for idx in 1:length(BETA_VALUES)]
survival = [results[idx]["survival"] / n_files[idx] for idx in 1:length(BETA_VALUES)]
total_extinctions = 1 .- survival

# Standard errors (using binomial proportion standard error formula for now)
se_demographic = sqrt.((demographic_extinctions .* (1 .- demographic_extinctions)) ./ n_files)
se_immune = sqrt.((immune_extinctions .* (1 .- immune_extinctions)) ./ n_files)
se_total = sqrt.((total_extinctions .* (1 .- total_extinctions)) ./ n_files)

# Plot the extinction probabilities
p = plot(BETA_VALUES, demographic_extinctions, ribbon=se_demographic, label="Demographic Extinction", alpha=0.5, linecolor=:blue)
plot!(p, BETA_VALUES, immune_extinctions, ribbon=se_immune, label="Immune Extinction", alpha=0.5, linecolor=:red)
plot!(p, BETA_VALUES, total_extinctions, ribbon=se_total, label="Total Extinction", alpha=0.5, linecolor=:green)
# plot!(p,BETA_VALUES,1.0 ./ BETA_VALUES .^ 0.2,linestyle=:dash)
xlabel!(p, "Beta Value")
ylabel!(p, "Probability of Extinction")
title!(p, "Extinction Probabilities vs. Beta")


# Save the plot
savefig("extinction_probabilities_vs_beta.png")
