using Plots
using Serialization
using Glob
using Statistics
using CSV
using DataFrames
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

# Constants
const DIRECTORY_PATH = "C:/Users/Daniel/Desktop/simresults_initial_variance"
const OUTPUT_PATH = "plotted_results_initial_variance/"
const variance_VALUES = LinRange(0.01, 0.1, 6)
const NUM_REPLICATES = 3000

# Ensure the output directory exists
isdir(OUTPUT_PATH) || mkdir(OUTPUT_PATH)

# Function to calculate the probability of survival
function calculate_survival_probability(variance)
    files = glob("simulation_results_variance_$(variance)_replicate_*.jld2", DIRECTORY_PATH)
    println(length(files))
    survival_count = 0
    total_variance_at_max = 0.0

    for file in files
        simulation = open(deserialize, file)
        total_infected = calculate_total_infected(simulation)
        antigenic_variance = calculate_antigenic_variance_per_deme(simulation)

        max_infected_index = argmax(total_infected)
        total_variance_at_max += antigenic_variance[1,max_infected_index]

        if total_infected[end] > 0
            survival_count += 1
        end
    end

    probability = survival_count / length(files)
    standard_error = sqrt(probability * (1 - probability) / length(files))
    return probability, standard_error, total_variance_at_max / length(files)
end

# Arrays to store probabilities and errors
probabilities = Float64[]
errors = Float64[]
avg_variances_at_max = Float64[]
# Calculate probabilities and errors for each variance
for variance in variance_VALUES
    prob, err, var = calculate_survival_probability(variance)
    push!(probabilities, prob)
    push!(errors, err)
    push!(avg_variances_at_max, var)
end

# Prepare the design matrix and response vector
X = hcat(ones(length(avg_variances_at_max)), avg_variances_at_max)
y = probabilities

# Solve the normal equations for linear regression
coefficients = X \ y
intercept, slope = coefficients[1], coefficients[2]
linear_fit = x -> slope * x + intercept

# Plotting
p = scatter(avg_variances_at_max, probabilities, ribbon=errors, label="Survival Probability", xlabel="Initial variance", ylabel="Probability", title="Survival Probability vs Initial variance", legend=:topright)
plot!(p, avg_variances_at_max, linear_fit.(avg_variances_at_max), label="Linear Fit", linestyle=:dash, linewidth=2, color=:red)
annotate!(p, [(maximum(avg_variances_at_max) * 0.5, maximum(probabilities) * 0.9, text("y = $(round(slope, digits=3))x + $(round(intercept, digits=3))", :left, 10, :red))])

# Save plot
savefig(p, "$(OUTPUT_PATH)/survival_probability_vs_variance.png")

# Print slope and intercept
println("Slope: ", slope)
println("Intercept: ", intercept)


# Save variances and probabilities to CSV
df = DataFrame(variance=avg_variances_at_max, probability=probabilities, error=errors)
CSV.write("../final_plots_viral_coev/variances_and_probabilities.csv", df)