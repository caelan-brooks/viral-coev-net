using Plots
using Serialization
using Glob
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

# Constants
const DIRECTORY_PATH = "C:/Users/Daniel/Desktop/simresults_initial_width"
const OUTPUT_PATH = "plotted_results_initial_width/"
const WIDTH_VALUES = LinRange(0.3, 6, 10)
const NUM_REPLICATES = 3000

# Ensure the output directory exists
isdir(OUTPUT_PATH) || mkdir(OUTPUT_PATH)

# Function to calculate the probability of survival
function calculate_survival_probability(width)
    files = glob("simulation_results_width_$(width)_replicate_*.jld2", DIRECTORY_PATH)
    println(length(files))
    survival_count = 0

    for file in files
        simulation = open(deserialize, file)
        total_infected = calculate_total_infected(simulation)
        if total_infected[end] > 0
            survival_count += 1
        end
    end

    probability = survival_count / length(files)
    standard_error = sqrt(probability * (1 - probability) / length(files))
    return probability, standard_error
end

# Arrays to store probabilities and errors
probabilities = []
errors = []

# Calculate probabilities and errors for each width
for width in WIDTH_VALUES
    prob, err = calculate_survival_probability(width)
    push!(probabilities, prob)
    push!(errors, err)
end

# Plotting
p = plot(WIDTH_VALUES, probabilities, ribbon=errors, label="Survival Probability", xlabel="Initial Width", ylabel="Probability", title="Survival Probability vs Initial Width", legend=:none)
savefig(p, "$(OUTPUT_PATH)/survival_probability_vs_width.png")
