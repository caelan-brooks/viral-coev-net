using Glob
using Serialization
using Statistics
using Plots
using Printf
using Base.Threads
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

# const DIRECTORY_PATH = "simresults_newseed/"

const DIRECTORY_PATH ="C:/Users/Daniel/Desktop/simresults/"
const OUTPUT_PATH = "plotted_results_new/"

const MIGRATION_RATES = exp10.(LinRange(-7.0, -0.5, 9)) # Example migration rates to sweep over

function plot_variance_vs_infected(migration_rate)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    println(length(files))
    short_migration_rate = @sprintf("%.2e", migration_rate)  # Format migration rate

    p = plot(;ylim=(0,1), xscale=:log10, xlabel="Total Infected", ylabel="Variance", title="MR: $short_migration_rate", legend=:none)

    # Iterate over files
    for file in files
        simulation = open(deserialize, file)
        # if calculate_total_infected(simulation)[end] == 0
        #     continue
        # end
        variances_per_deme = calculate_antigenic_variance_per_deme(simulation)
        total_infected_per_deme = calculate_total_infected_per_deme(simulation)

        num_demes = size(variances_per_deme, 1)

        # Iterate over demes
        for deme in 1:num_demes
            max_index = findfirst(==(maximum(total_infected_per_deme[deme, :])), total_infected_per_deme[deme, :])

            # Create the Boolean mask
            reg = (total_infected_per_deme[deme, :] .> 0) .& 
                (variances_per_deme[deme, :] .> 0) .& 
                (1:length(total_infected_per_deme[deme, :]) .< max_index)

            # Filter data using the mask
            if sum(reg) > 0
                filtered_infected = total_infected_per_deme[deme, reg]
                filtered_variances = variances_per_deme[deme, reg]

                # Plotting with see-through markers
                scatter!(p, filtered_infected, filtered_variances, m=(:cross, 2, 0.5), color=(deme == 1 ? :blue : :red))
            end
        end
    end

    # Save the plot
    # savefig(p, "$(OUTPUT_PATH)/variance_vs_infected_mr$(idx).png")
    return p
end

# Main loop to create a 3x3 grid of plots
all_plots = []
for (idx, migration_rate) in enumerate(MIGRATION_RATES)
    println("Processing migration rate: $migration_rate")
    p = plot_variance_vs_infected(migration_rate)
    push!(all_plots, p)
end

# Combine all plots into a 3x3 grid layout
combined_plot = plot(all_plots..., layout=(3, 3), size=(900, 900))

# Save the combined plot
savefig(combined_plot, "$(OUTPUT_PATH)/combined_variance_vs_infected.png")