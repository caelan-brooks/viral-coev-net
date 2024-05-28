using Base.Threads
using Serialization
using Random
using LinearAlgebra
using CSV
using DataFrames
include("/home/dswartz/viral-coev-net/julia_code/coevolution_network_base.jl")
using .CoevolutionNetworkBase

const OUTPUT_DIRECTORY = "/pool001/dswartz/real_network_PL_with_dx"

println("Number of threads: ", nthreads())

const N0 = 100
const L = 50.0  # Consider if this needs to be longer
const dx = 0.1
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const beta = 2.5
const alpha = 0.0
const gamma = 1.0
const D = 0.0042 # make this smaller? 
const sigma = 2.0
const DURATION = 80.0
const DT = 0.05
const THIN_BY = 20
const NUM_REPLICATES = 4000
const START_REPLICATE = 1

df = CSV.read("cleaned_adjacency_matrix.csv", DataFrame)
migration_matrix = Matrix(df[:, 1:end-1])  # Assuming the last column is the population sizes vector
population_sizes = df[:,:population_sizes]

println(migration_matrix)
println(population_sizes)

function run_single_simulation(args)
    # Unpack arguments
    outbreak_deme_idx, simulation_number = args
    println(simulation_number)
    flush(stdout)

    # Generate unique seed and set random seed
    seed = hash((outbreak_deme_idx, simulation_number))
    Random.seed!(seed)

    network_size = length(population_sizes)

    # Initialize viral and immune densities
    viral_densities = [zeros(Float64, length(x)) for _ in 1:network_size]
    immune_densities = [zeros(Float64, length(x)) for _ in 1:network_size]

    # Set the value of viral_density at the closest index to N0/dx for the first population
    index_closest_to_zero = argmin(abs.(x))
    viral_densities[outbreak_deme_idx][index_closest_to_zero] = N0/dx

    # Create populations
    
    populations = [Population(L, dx, r, M, beta, alpha, gamma, D, population_sizes[i], viral_densities[i], immune_densities[i]; sigma=sigma, noise_method=:PL_with_dx) for i in 1:network_size]

    # Initialize populations and network with the new migration matrix
    network = Network(populations, migration_matrix)
    simulation = Simulation(network, DT, DURATION; thin_by=THIN_BY)

   # Run the simulation
    try
        @time run_simulation!(simulation)

        # Calculate total infected per deme
        total_infected_per_deme = calculate_total_infected_per_deme(simulation)
        antigenic_variance_per_deme = calculate_antigenic_variance_per_deme(simulation)

        # Prepare output file path using migration rate index
        output_file = joinpath(OUTPUT_DIRECTORY, "outbreak_deme_idx_$(outbreak_deme_idx)", "replicate_$(simulation_number).jld2")

        # Save the calculated total infected per deme
        open(output_file, "w") do file
            serialize(file, (total_infected_per_deme, antigenic_variance_per_deme))
        end
    catch e
        println("Error in simulation with args $(args): $e")
    end

end

# Function to save parameters and migration rates to CSV
function save_parameters_and_migration_rates_to_csv(output_directory)
    # Create a DataFrame for the main parameters
    parameters = DataFrame(
        L = [L],
        dx = [dx],
        r = [r],
        M = [M],
        beta = [beta],
        alpha = [alpha],
        gamma = [gamma],
        D = [D],
        DURATION = [DURATION],
        DT = [DT],
        THIN_BY = [THIN_BY],
        NUM_REPLICATES = [NUM_REPLICATES],
        START_REPLICATE = [START_REPLICATE]
    )
    
    # Save the main parameters DataFrame to a CSV file
    parameters_path = joinpath(output_directory, "simulation_parameters.csv")
    CSV.write(parameters_path, parameters)
end

save_parameters_and_migration_rates_to_csv(OUTPUT_DIRECTORY)

function main(job_id_arg)
    job_id = parse(Int, job_id_arg)
    total_combinations = length(population_sizes)
    
    # Check if job_id is within the valid range
    if job_id < 1 || job_id > total_combinations
        error("Invalid job ID: $job_id")
    end

    outbreak_deme_idx = job_id

    # Prepare output directory for this combination of parameters
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "outbreak_deme_idx_$(outbreak_deme_idx)")
    mkpath(output_subdirectory)

    # Run simulations for all replicates
    @threads for simulation_number in START_REPLICATE:(START_REPLICATE+NUM_REPLICATES)
        args = (outbreak_deme_idx, simulation_number)
        run_single_simulation(args)
    end
end

# Get job_id from command line arguments or environment variable
job_id_arg = get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
main(job_id_arg)