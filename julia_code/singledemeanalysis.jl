using Serialization
using Glob
using Statistics
using LaTeXStrings
using Plots
const DIRECTORY_PATH = "simresults_singledeme"

function calculate_probability_of_survival(cutoff::Real)
    files = glob("simulation_results_replicate_*.jld2", DIRECTORY_PATH)
    
    survival_counts = 0
    total_files = length(files)
    println("Total files: ", total_files)

    plot(;  xlabel="Time", 
            ylabel="Number of Infected Individuals", 
            title="Infection Over Time",
            yscale=:log10,
            guidefont=font(16, "Computer Modern"), 
            tickfont=font(14, "Computer Modern"), 
            legendfont=font(14, "Computer Modern"), 
            titlefont=font(16, "Computer Modern")
        )

    for file in files
        result_dict = open(deserialize, file)
        xdata = result_dict["times"]
        ydata = result_dict["total_infected_number"]
        reg = ydata .> 0.0
        plot!(xdata[reg], ydata[reg], legend=false)

        data_slice = result_dict["total_infected_number"][Int(round(end * 0.95)) : end]
        
        if mean(data_slice) > cutoff
            survival_counts += 1
        end
    end

    
    probability_of_survival = survival_counts / total_files
    standard_error = sqrt((probability_of_survival * (1 - probability_of_survival)) / total_files)
    
    return probability_of_survival, standard_error
end

# Usage:
# Adjust the cutoff value as per your requirement
cutoff = 0  # Adjust this value

prob_of_survival, std_error = calculate_probability_of_survival(cutoff)
println("Probability of Survival: ", prob_of_survival)
println("Standard Error: ", std_error)
display(current())
