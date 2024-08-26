using Base.Threads
using Serialization
using Random
using LinearAlgebra
using CSV
using DataFrames
include("/home/dswartz/viral-coev-net/julia_code/coevolution_network_base.jl")
using .CoevolutionNetworkBase


const OUTPUT_DIRECTORY = "/pool001/dswartz/peak_scaling"
const MIGRATION_RATES = exp10.(LinRange(-8.0, -2, 12))
const HOST_POPULATION_SIZES = exp10.(4:0.5:8)

println("Number of threads: ", nthreads())

const N0 = 100
const L = 40.0
const dx = 0.1
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const beta = 2.5
const alpha = 0.0
const gamma = 1.0
const D = 0.01
const sigma = 2.0 
const DURATION = 100.0
const DT = 0.05
const THIN_BY = 10
const NUM_REPLICATES = 10000
const START_REPLICATE = 1 
const noise_method = :PL_with_dx

function run_single_simulation(args)
    # Unpack arguments
    host_per_deme_idx, migration_rate_idx, simulation_number = args
    println(simulation_number)
    flush(stdout)
    
    # Retrieve the actual migration rate using the index
    migration_rate = MIGRATION_RATES[migration_rate_idx]
    host_population_per_deme = HOST_POPULATION_SIZES[host_per_deme_idx]

    # Generate unique seed and set random seed
    seed = hash((host_per_deme_idx, migration_rate_idx, simulation_number))
    Random.seed!(seed)

    # Retrieve the adjacency matrix for the given index
    adjacency_matrix = [0.0 1.0; 1.0 0.0]
    network_size = 2

    # Initialize viral and immune densities
    viral_densities = [zeros(Float64, length(x)) for _ in 1:network_size]
    immune_densities = [zeros(Float64, length(x)) for _ in 1:network_size]

    # Set the value of viral_density at the closest index to N0/dx for the first population
    index_closest_to_zero = argmin(abs.(x))
    viral_densities[1][index_closest_to_zero] = N0/dx
    # initial_antigenic_variance = 0.1;
    # viral_densities[1] .= N0/sqrt(2 * pi * initial_antigenic_variance) .* exp.(-x.^2/2/initial_antigenic_variance)

    # Create populations
    populations = [Population(L, dx, r, M, beta, alpha, gamma, D, host_population_per_deme, viral_densities[i], immune_densities[i]; sigma=sigma, noise_method=noise_method) for i in 1:network_size]

    # Initialize populations and network with the new migration matrix
    migration_matrix = migration_rate * adjacency_matrix
    network = Network(populations, migration_matrix)
    simulation = Simulation(network, DT, DURATION; thin_by=THIN_BY)

   # Run the simulation
    try
        @time run_simulation!(simulation);
        # thin_simulation!(simulation, THIN_BY)

        # Calculate total infected per deme
        total_infected_per_deme = calculate_total_infected_per_deme(simulation)
        antigenic_variance_per_deme = calculate_antigenic_variance_per_deme(simulation)

        # Prepare output file path using migration rate index
        output_file = joinpath(OUTPUT_DIRECTORY, "migration_rate_idx_$(migration_rate_idx)", "replicate_$(simulation_number).jld2")

        # Save the calculated total infected per deme
        open(output_file, "w") do file
            serialize(file, (total_infected_per_deme, antigenic_variance_per_deme, simulation.duration_times))
        end
    catch e
        println("Error in simulation with args $(args): $e")
    end

end

# Function to save parameters and migration rates to CSV
function save_parameters_and_migration_rates_to_csv(output_directory, migration_rates, host_population_sizes)
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

    population_sizes_df = DataFrame(HOST_POPULATION_SIZES = host_population_sizes)
    population_sizes_path = joinpath(output_directory, "host_population_sizes.csv")
    CSV.write(population_sizes_path, population_sizes_df)
end

save_parameters_and_migration_rates_to_csv(OUTPUT_DIRECTORY, MIGRATION_RATES)

function main(job_id_arg)
    job_id = parse(Int, job_id_arg)
    total_combinations = length(MIGRATION_RATES) * length(HOST_POPULATION_SIZES)
    
    # Check if job_id is within the valid range
    if job_id < 1 || job_id > total_combinations
        error("Invalid job ID: $job_id")
    end

    host_per_deme_idx = (job_id - 1) ÷ length(MIGRATION_RATES) + 1
    migration_rate_idx = (job_id - 1) % length(MIGRATION_RATES) + 1

    # Prepare output directory for this combination of parameters
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "host_per_deme_idx_$(host_per_deme_idx)", "migration_rate_idx_$(migration_rate_idx)")
    mkpath(output_subdirectory)

    # Run simulations for all replicates
    @threads for simulation_number in START_REPLICATE:(START_REPLICATE+NUM_REPLICATES-1)
        args = (host_per_deme_idx, migration_rate_idx, simulation_number)
        run_single_simulation(args)
        GC.gc() # Careful, this might not be the right place to clear memory
    end
end

# Get job_id from command line arguments or environment variable
job_id_arg = get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
main(job_id_arg)