using Glob
using Serialization
using Plots
using Printf
using Statistics

include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DIRECTORY_PATH = "C:/Users/Daniel/Desktop/simresults_random_rates/"
const OUTPUT_PATH = "plotted_results_random_rates/"

function analyze_and_plot(idx, migration_rate)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    println(length(files))

    # Create empty vectors to hold the densities for each deme and type
    viral_densities_deme1 = Vector{Float64}[]
    immune_densities_deme1 = Vector{Float64}[]
    viral_densities_deme2 = Vector{Float64}[]
    immune_densities_deme2 = Vector{Float64}[]

    simulation = open(deserialize,files[1])
    xs = simulation.trajectory[end].populations[1].xs
    Nh = simulation.trajectory[end].populations[1].Nh
    
    Threads.@threads for file in files
        simulation = open(deserialize, file)
        
        # Check if total viral density is non-zero at the end of the simulation for deme 1 and 2
        end_state = simulation.trajectory[end]
        
        total_viral_deme1 = sum(end_state.populations[1].viral_density .* end_state.populations[1].dx)
        # println(total_viral_deme1)
        total_viral_deme2 = sum(end_state.populations[2].viral_density .* end_state.populations[2].dx)

        if total_viral_deme1 > 0.0
            push!(viral_densities_deme1, end_state.populations[1].viral_density)
            push!(immune_densities_deme1, end_state.populations[1].immune_density)
        end

        if total_viral_deme2 > 0.0
            push!(viral_densities_deme2, end_state.populations[2].viral_density)
            push!(immune_densities_deme2, end_state.populations[2].immune_density)
        end
        
    end

    p = plot(layout=(2, 2), size=(800, 600), title="Migration Rate: $migration_rate")

    # Plot densities for Deme 1
    for vd in viral_densities_deme1
        plot!(p[1, 1], xs, vd, label=false, legend=false, color=:blue, alpha=0.5, ylim=(0, Nh))
    end
    # plot!(p[1, 1], xlabel="Antigenic Coordinate \$x\$", title="Viral Density in Deme 1")
    
    for id in immune_densities_deme1
        plot!(p[2, 1], xs, id, label=false, legend=false, color=:red, alpha=0.5)
    end
    plot!(p[2, 1], xlabel="Antigenic Coordinate \$x\$", title="Immune Density in Deme 1")
    
    # Plot densities for Deme 2
    for vd in viral_densities_deme2
        plot!(p[1, 2], xs, vd, label=false, legend=false, color=:blue, alpha=0.5, ylim=(0, Nh))
    end
    plot!(p[1, 2], xlabel="Antigenic Coordinate \$x\$", title="Viral Density in Deme 2")
    
    for id in immune_densities_deme2
        plot!(p[2, 2], xs, id, label=false, legend=false, color=:red, alpha=0.5)
    end
    plot!(p[2, 2], xlabel="Antigenic Coordinate \$x\$", title="Immune Density in Deme 2")
    

    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    savefig(p, "$(OUTPUT_PATH)/densities_idx_$(idx)_migration_$(migration_rate).png")
end

function main()
    # migration_rates = vcat([0], exp10.(LinRange(-6, 0.5, 9))) # Example migration rates to sweep over
    migration_rates = exp10.(LinRange(-6,0.5,9))
    for (idx, migration_rate) in enumerate(migration_rates)
        println("Analyzing migration rate: $migration_rate")
        analyze_and_plot(idx,migration_rate)
    end
end

main()
