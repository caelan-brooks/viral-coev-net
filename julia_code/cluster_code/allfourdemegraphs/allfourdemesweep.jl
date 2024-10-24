using Base.Threads
using Serialization
using Random
using LinearAlgebra
using DataFrames
using CSV
include("/home/dswartz/viral-coev-net/julia_code/coevolution_network_base.jl")
using .CoevolutionNetworkBase

# New line to include the script for loading adjacency matrices
include("./all_adjacency_matrices.jl")

const OUTPUT_DIRECTORY = "/pool001/dswartz/all_four_deme_PL_with_dx"
const MIGRATION_RATES = exp10.(LinRange(-8, 0.0, 10))

println("Number of threads: ", nthreads())

const TOTAL_HOST_POPULATION = 4 * 10^6
const L = 40.0
const dx = 0.1
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const beta = 2.5
const alpha = 1.0
const gamma = 0.0
const D = 0.01
const sigma = 2.0
const noise_method = :PL_with_dx
const DURATION = 80.0
const DT = 0.05
const THIN_BY = 20
const NUM_REPLICATES = 9000
const START_REPLICATE = 1001

function run_single_simulation(args)
    # Unpack arguments
    migration_rate_idx, adjacency_matrix_idx, simulation_number = args
    println(simulation_number)
    flush(stdout)
    
    # Retrieve the actual migration rate using the index
    migration_rate = MIGRATION_RATES[migration_rate_idx]

    # Generate unique seed and set random seed
    seed = hash((migration_rate_idx, adjacency_matrix_idx, simulation_number))
    Random.seed!(seed)

    # Retrieve the adjacency matrix for the given index
    adjacency_matrix = all_adjacency_matrices[adjacency_matrix_idx]
    network_size = size(adjacency_matrix, 1)

    # Calculate population per deme
    population_per_deme = round(Int, TOTAL_HOST_POPULATION / network_size)

    # Initialize viral and immune densities
    viral_densities = [zeros(Float64, length(x)) for _ in 1:network_size]
    immune_densities = [zeros(Float64, length(x)) for _ in 1:network_size]

    # Set the value of viral_density at the closest index to 1/dx for the first population
    index_closest_to_zero = argmin(abs.(x))
    viral_densities[1][index_closest_to_zero] = 100/dx

    # Create populations
    populations = [Population(L, dx, r, M, beta, alpha, gamma, D, population_per_deme, viral_densities[i], immune_densities[i]; sigma = sigma, noise_method=noise_method) for i in 1:network_size]

    # Initialize populations and network with the new migration matrix
    migration_matrix = migration_rate * adjacency_matrix
    network = Network(populations, migration_matrix)
    simulation = Simulation(network, DT, DURATION; thin_by=THIN_BY)

   # Run the simulation
    try
        @time run_simulation!(simulation)
        # thin_simulation!(simulation, THIN_BY)

        # Calculate total infected per deme
        total_infected_per_deme = calculate_total_infected_per_deme(simulation)
        # antigenic_variance_per_deme = calculate_antigenic_variance_per_deme(simulation)

        # Prepare output file path using migration rate index
        output_file = joinpath(OUTPUT_DIRECTORY, "adjacency_matrix_idx_$(adjacency_matrix_idx)", "migration_rate_idx_$(migration_rate_idx)", "replicate_$(simulation_number).jld2")

        # Save the calculated total infected per deme
        open(output_file, "w") do file
            # serialize(file, (total_infected_per_deme, antigenic_variance_per_deme))
            serialize(file, total_infected_per_deme)
        end
    catch e
        println("Error in simulation with args $(args): $e")
    end

end

# Function to save parameters and migration rates to CSV
function save_parameters_and_migration_rates_to_csv(output_directory, migration_rates)
    # Create a DataFrame for the main parameters
    parameters = DataFrame(
        TOTAL_HOST_POPULATION = [TOTAL_HOST_POPULATION],
        L = [L],
        dx = [dx],
        r = [r],
        M = [M],
        beta = [beta],
        alpha = [alpha],
        gamma = [gamma],
        D = [D],
        sigma = [sigma],
        DURATION = [DURATION],
        DT = [DT],
        THIN_BY = [THIN_BY],
        NUM_REPLICATES = [NUM_REPLICATES],
        START_REPLICATE = [START_REPLICATE],
        noise_method = [noise_method]
    )
    
    # Save the main parameters DataFrame to a CSV file
    parameters_path = joinpath(output_directory, "simulation_parameters.csv")
    mkpath(output_directory)
    CSV.write(parameters_path, parameters)
    
    # Create a DataFrame for migration rates
    migration_rates_df = DataFrame(MIGRATION_RATES = migration_rates)
    
    # Save the migration rates DataFrame to a CSV file
    migration_rates_path = joinpath(output_directory, "migration_rates.csv")
    CSV.write(migration_rates_path, migration_rates_df)
end

save_parameters_and_migration_rates_to_csv(OUTPUT_DIRECTORY, MIGRATION_RATES)

function main(job_id_arg)
    job_id = parse(Int, job_id_arg)
    total_combinations = length(all_adjacency_matrices) * length(MIGRATION_RATES)
    
    # Check if job_id is within the valid range
    if job_id < 1 || job_id > total_combinations
        error("Invalid job ID: $job_id")
    end

    # Calculate indices for adjacency matrix and migration rate
    adjacency_matrix_idx = (job_id - 1) ÷ length(MIGRATION_RATES) + 1
    migration_rate_idx = (job_id - 1) % length(MIGRATION_RATES) + 1

    # Prepare output directory for this combination of parameters
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "adjacency_matrix_idx_$(adjacency_matrix_idx)", "migration_rate_idx_$(migration_rate_idx)")
    mkpath(output_subdirectory)

    # Run simulations for all replicates
    @threads for simulation_number in START_REPLICATE:(START_REPLICATE+NUM_REPLICATES-1)
        args = (migration_rate_idx, adjacency_matrix_idx, simulation_number)
        run_single_simulation(args)
    end
end

# Get job_id from command line arguments or environment variable
job_id_arg = get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
main(job_id_arg)

