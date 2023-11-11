using Glob
using Plots
using Serialization
using Statistics

# Constants and parameters setup
cutoff = 100
migration_rates = exp10.(LinRange(-5, 0, 10))
mutation_rates = exp10.(LinRange(-6, -5, 3))
results_dir = "simulation_results"

# Function to calculate the probability for strain two exceeding the cutoff
function calculate_probability_for_file(file_path)
    trajectories = deserialize(file_path)
    exceed_count = sum(any(trajectory[:, 2, :] .> cutoff) for trajectory in trajectories)
    exceed_probability = exceed_count / length(trajectories)
    return exceed_probability
end



# Loop through all combinations of rate indices
for (mut_idx, μ_rate) in enumerate(mutation_rates)
    probabilities = []

    # Plotting setup
    plt = plot(legend=:bottomleft)  # Initialize an empty plot
    xlabel!("Migration Rate k")
    ylabel!("Probability of Strain 2 Exceeding $cutoff")
    title!("Probability vs Migration Rate")
    for (k_idx, k_rate) in enumerate(migration_rates)
        # Generate the expected filename based on current rate indices
        file_path = joinpath(results_dir, "trajectories_k$(k_idx)_μ$(mut_idx).dat")
        if isfile(file_path)
            println("yay")
            prob = calculate_probability_for_file(file_path)
            push!(probabilities, prob)
        else
            println("nay")
            push!(probabilities, NaN)  # If file does not exist, use NaN to denote missing data
        end
    end
    println(probabilities)
    # Add the current mutation rate's probabilities to the plot
    plot!(plt, migration_rates, probabilities, label="μ = $(μ_rate)", xscale=:log10)
    # Save the figure
    savefig(plt,"Probability_vs_Migration_Rate_$(μ_rate).png")
end


