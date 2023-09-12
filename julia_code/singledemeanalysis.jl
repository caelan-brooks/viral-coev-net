using Serialization
using Glob

const DIRECTORY_PATH = "simresults_singledeme"

function calculate_probability_of_survival(cutoff)
    files = glob("simulation_results_replicate_*.jld2", DIRECTORY_PATH)
    
    survival_counts = 0
    total_files = length(files)
    println("Total files: ", total_files)
    
    for file in files
        result_dict = open(deserialize, file)
        
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
cutoff = 10^2  # Adjust this value

prob_of_survival, std_error = calculate_probability_of_survival(cutoff)
println("Probability of Survival: ", prob_of_survival)
println("Standard Error: ", std_error)
