using Glob
using Serialization
using Statistics
using Plots

const DIRECTORY_PATH = "simresults_largedt"

function calculate_probability_of_survival(migration_rate, cutoff)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    survival_counts = 0
    total_files = length(files)
    println(total_files)
    
    for file in files
        result_dict = open(deserialize, file) # Changed this line
        
        data_slice = result_dict["total_infected_number"][Int(round(end * 0.95)) : end]
        
        if mean(data_slice) > cutoff
            survival_counts += 1
        end
    end
    
    probability_of_survival = survival_counts / total_files
    standard_error = sqrt((probability_of_survival * (1 - probability_of_survival)) / total_files)
    
    return probability_of_survival, standard_error
end

function plot_total_infected_trajectories(migration_rate)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    plot()
    
    for file in files
        result_dict = open(deserialize, file) # Changed this line
        
        plot!(result_dict["times"], result_dict["total_infected_number"], label="Replicate", yscale=:log10)
    end
    
    xlabel!("Time")
    ylabel!("Total Infected Number")
    title!("Total Infected Number vs Time (Migration Rate: $migration_rate)")
    display(plot())
end


function main()
    migration_rates = exp10.(LinRange(-6, 0.5, 9))
    cutoff = 10^3
    
    probabilities = []
    errors = []
    
    for migration_rate in migration_rates
        prob, err = calculate_probability_of_survival(migration_rate, cutoff)
        push!(probabilities, prob)
        push!(errors, err)
        
        # Uncomment the following line to plot the trajectories for each migration rate
        # plot_total_infected_trajectories(migration_rate)
    end

    plotvar = plot(migration_rates, probabilities, yerr=errors, seriestype=:scatter, xscale=:log10,
         xlabel="Migration Rate", ylabel="Probability of Survival",
         title="Probability of Survival as a Function of Migration Rate",
         legend=false)
         
    display(plotvar)
    
    # Save the plot as a PNG file
    savefig(plotvar, "probability_of_survival.png")
end
main()
