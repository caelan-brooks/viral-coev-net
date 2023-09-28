using Glob
using Serialization
using Statistics
using Plots
using Printf
using Base.Threads
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

# const DIRECTORY_PATH = "simresults_newseed/"

const DIRECTORY_PATH ="C:/Users/Daniel/Desktop/simresults_oneinfected"
const OUTPUT_PATH = "plotted_results_oneinfected"

function calculate_total_infected_per_deme(simulation::Simulation)
    # Get the number of populations (demes)
    num_demes = length(simulation.trajectory[1].populations)
    
    # Initialize a list of vectors to hold the total number of infected individuals at each time step for each deme
    total_infected_per_deme = [zeros(Float64, length(simulation.trajectory)) for _ in 1:num_demes]

    # Iterate over each network state in the simulation's trajectory
    for (i, network) in enumerate(simulation.trajectory)
        # Iterate over each population in the current network state
        for (j, population) in enumerate(network.populations)
            # Increment the total infected count for the current time step and deme
            # by adding the integral of the viral density (viral_density multiplied by dx)
            total_infected_per_deme[j][i] += sum(population.viral_density .* population.dx)
        end
    end

    # Return the list containing the total number of infected individuals at each time step for each deme
    return total_infected_per_deme
end

function plot_total_infected_per_deme(migration_rate, migration_rate_index)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    formatted_migration_rate = @sprintf("%.2e", migration_rate)

    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    plots_surv = [plot(title="Surviving Trajectories (Migration Rate: $formatted_migration_rate, Deme: $i)", xlabel="Time", ylabel="Total Infected Number in Deme $i", legend=false, ylims=(1,3*10^6)) for i in 1:2]
    plots_not_surv = [plot(title="Non-Surviving Trajectories (Migration Rate: $formatted_migration_rate, Deme: $i)", xlabel="Time", ylabel="Total Infected Number in Deme $i", legend=false, ylims=(1,3*10^6)) for i in 1:2]

    for file in files
        simulation = open(deserialize, file)
        xdata = simulation.duration_times
        ydata_list = calculate_total_infected_per_deme(simulation) # Updated function call

        for (idx, ydata) in enumerate(ydata_list)
            reg = ydata .> 0.0

            if ydata[end] > 0.0
                file_suffix = "surv"
                if any(reg)
                    plot!(plots_surv[idx], xdata[reg], ydata[reg], label=false, legend=false, yscale=:log10)
                else
                    plot!(plots_surv[idx], xdata[reg], ydata[reg], label=false, legend=false)
                end
            else
                file_suffix = "not_surv"
                if any(reg)
                    plot!(plots_not_surv[idx], xdata[reg], ydata[reg], label=false, legend=false, yscale=:log10)
                else
                    plot!(plots_not_surv[idx], xdata[reg], ydata[reg], label=false, legend=false)
                end
            end
        end
    end

    for (idx, p) in enumerate(plots_surv)
        savefig(p, "$(OUTPUT_PATH)/trajectory_surv_deme_$(idx)_idx_$(migration_rate_index)_rate_$(formatted_migration_rate).png")
    end
    for (idx, p) in enumerate(plots_not_surv)
        savefig(p, "$(OUTPUT_PATH)/trajectory_not_surv_deme_$(idx)_idx_$(migration_rate_index)_rate_$(formatted_migration_rate).png")
    end
end




function calculate_probability_of_survival(migration_rate, cutoff)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    survival_counts = 0
    total_files = length(files)
    println(total_files)
    
    for file in files
            simulation = open(deserialize, file) # Changed this line
            
            data_slice = calculate_total_infected(simulation)[Int(round(end * 0.95)) : end]
            
            if mean(data_slice) > cutoff
                survival_counts += 1
            end
        end
    
    probability_of_survival = survival_counts / total_files
    standard_error = sqrt((probability_of_survival * (1 - probability_of_survival)) / total_files)
    
    return probability_of_survival, standard_error
end

function plot_total_infected_trajectories(migration_rate, migration_rate_index)
    files = glob("simulation_results_migration_$(migration_rate)_replicate_*.jld2", DIRECTORY_PATH)
    
    formatted_migration_rate = @sprintf("%.2e", migration_rate)

    plot_surviving = plot(title = "Surviving Trajectories (Migration Rate: $formatted_migration_rate)")
    plot_not_surviving = plot(title = "Non-Surviving Trajectories (Migration Rate: $formatted_migration_rate)")
    
    for file in files
        simulation = open(deserialize, file)
        xdata = simulation.duration_times
        ydata = calculate_total_infected(simulation)
        reg = ydata .> 0.0

        xlabel!("Time")
        ylabel!("Total Infected Number")

        if ydata[end] > 0.0
            plot!(plot_surviving, xdata[reg], ydata[reg], label=false, legend=false, yscale=:log10)
        else
            plot!(plot_not_surviving, xdata[reg], ydata[reg], label=false, legend=false, yscale=:log10)
        end
    end

    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    savefig(plot_surviving, "$(OUTPUT_PATH)/trajectory_surv_idx_$(migration_rate_index)_rate_$(formatted_migration_rate).png")
    savefig(plot_not_surviving, "$(OUTPUT_PATH)/trajectory_not_surv_idx_$(migration_rate_index)_rate_$(formatted_migration_rate).png")
end



function main()
    println("Number of threads: ", Threads.nthreads())

    migration_rates = vcat([0], exp10.(LinRange(-6, 0.5, 9))) # Example migration rates to sweep over
    cutoff = 10

    num_rates = length(migration_rates)
    probabilities = Vector{Float64}(undef, num_rates)
    errors = Vector{Float64}(undef, num_rates)

    # Use the @threads macro to run the loop in parallel
    for idx in 1:num_rates
        migration_rate = migration_rates[idx]
        prob, err = calculate_probability_of_survival(migration_rate, cutoff)

        # Assign the results to the specific index
        probabilities[idx] = prob
        errors[idx] = err

        # Uncomment the following lines to plot the trajectories for each migration rate
        plot_total_infected_trajectories(migration_rate, idx)
        plot_total_infected_per_deme(migration_rate, idx)
    end

    plotvar = plot(migration_rates[2:end], probabilities[2:end], yerr=errors[2:end], seriestype=:scatter, xscale=:log10,
         xlabel="Migration Rate", ylabel="Probability of Survival",
         title="Probability of Survival as a Function of Migration Rate",
         legend=false)
    
    # hline!([0.4185], color=:red, linestyle=:dash)
    # hline!([probabilities[1]], color=:blue, linestyle=:dash)
    display(plotvar)

    # Check if the "plotted_results" folder exists, if not create it
    if !isdir("$(OUTPUT_PATH)")
        mkdir("$(OUTPUT_PATH)")
    end

    # Save the plot as a PNG file in the "plotted_results" folder
    savefig(plotvar, "$(OUTPUT_PATH)/probability_of_survival.png")
end

main()
