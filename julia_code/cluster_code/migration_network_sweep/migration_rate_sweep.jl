using Base.Threads
using Serialization
using Random
using LinearAlgebra
include("/home/dswartz/viral-coev-net/julia_code/coevolution_network_base.jl")
using .CoevolutionNetworkBase

const OUTPUT_DIRECTORY = "/pool001/dswartz/migration_network_sweep_results"
const MIGRATION_RATES = exp10.(LinRange(-6, 0.5, 10))
const NETWORK_SIZES = 2:6

println("Number of threads: ", nthreads())

const TOTAL_HOST_POPULATION = 12 * 10^6
const L = 40.0
const dx = 0.3
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const beta = 2.5
const alpha = 1.0
const gamma = 0.0
const D = 0.01
const DURATION = 80.0
const DT = 0.05
const THIN_BY = 50
const NUM_REPLICATES = 3000

function run_single_simulation(args)
    # Unpack arguments
    migration_rate_idx, network_size, simulation_number = args
    println(simulation_number)
    flush(stdout)
    
    # Retrieve the actual migration rate using the index
    migration_rate = MIGRATION_RATES[migration_rate_idx]

    # Generate unique seed and set random seed
    seed = hash((migration_rate_idx, network_size, simulation_number))
    Random.seed!(seed)
    
    # Calculate normalized migration rate and population per deme
    normalized_migration_rate = migration_rate / (network_size - 1)
    population_per_deme = round(Int, TOTAL_HOST_POPULATION / network_size)

    
    # Initialize viral and immune densities
    viral_densities = [zeros(Float64, length(x)) for _ in 1:network_size]
    immune_densities = [zeros(Float64, length(x)) for _ in 1:network_size]

    # Set the value of viral_density at the closest index to 1/dx for the first population
    index_closest_to_zero = argmin(abs.(x))
    viral_densities[1][index_closest_to_zero] = 100/dx

    # Create populations
    populations = [Population(L, dx, r, M, beta, alpha, gamma, D, population_per_deme, viral_densities[i], immune_densities[i]) for i in 1:network_size]

    # Initialize populations and network
    migration_matrix = normalized_migration_rate * (ones(network_size, network_size) - I)
    network = Network(populations, migration_matrix)
    simulation = Simulation(network, DT, DURATION)

    # Run the simulation
    try
        @time run_simulation!(simulation)
        thin_simulation!(simulation, THIN_BY)

        # Prepare output file path using migration rate index
        output_file = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(migration_rate_idx)", "replicate_$(simulation_number).jld2")

        # Save the simulation results
        open(output_file, "w") do file
            serialize(file, simulation)
        end
    catch e
        println("Error in simulation with args $(args): $e")
    end
end

function main(job_id_arg)
    job_id = parse(Int, job_id_arg)
    total_combinations = length(NETWORK_SIZES) * length(MIGRATION_RATES)
    
    # Check if job_id is within the valid range
    if job_id < 1 || job_id > total_combinations
        error("Invalid job ID: $job_id")
    end

    # Calculate indices for network size and migration rate
    network_size_idx = (job_id - 1) ÷ length(MIGRATION_RATES) + 1
    migration_rate_idx = (job_id - 1) % length(MIGRATION_RATES) + 1
    network_size = NETWORK_SIZES[network_size_idx]

    # Prepare output directory for this combination of parameters
    output_subdirectory = joinpath(OUTPUT_DIRECTORY, "network_size_$(network_size)", "migration_rate_idx_$(migration_rate_idx)")
    mkpath(output_subdirectory)

    # Run simulations for all replicates
    @threads for simulation_number in 1:NUM_REPLICATES
        args = (migration_rate_idx, network_size, simulation_number)
        run_single_simulation(args)
    end
end



# Get job_id from command line arguments or environment variable
job_id_arg = get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
main(job_id_arg)

