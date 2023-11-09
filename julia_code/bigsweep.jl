using Distributions
using Base.Threads
using Serialization
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const OUTPUT_DIRECTORY = "C:/Users/Daniel/Desktop/simresults_random_rates"

if !isdir(OUTPUT_DIRECTORY)
    mkdir(OUTPUT_DIRECTORY)
end

println("Number of threads: ", nthreads())

const L = 40.0
const dx = 0.3
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const beta = 2.5
const alpha = 1.0
const gamma = 0.0
const D = 0.01
# const D = 0.0001
const Nh = 3 * 10^6
const dt = 0.05
const duration = 80.0

# Find the index of the grid point closest to x = 0
index_closest_to_zero = argmin(abs.(x))

# Initialize viral_density as zeros
const viral_density = zeros(Float64, length(x))

# Set the value of viral_density at the closest index to 1/dx
viral_density[index_closest_to_zero] = 1/dx

const viral_density2 = zeros(Float64, length(x))
const immune_density = zeros(Float64, length(x))

const thin_by = 50

using Random

function run_single_simulation(args)
    migration_rate, simulation_number = args

    # Generating a unique seed using migration_rate and simulation_number
    seed = hash((migration_rate, simulation_number))
    Random.seed!(seed)

    population1 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density)
    population2 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density2, immune_density)
    
    # migration_matrix = [0.0 migration_rate; migration_rate 0.0]
    # Generate off-diagonal terms as iid exponential random variables
    migration1to2 = rand(Exponential(migration_rate))
    migration2to1 = rand(Exponential(migration_rate))
    migration_matrix = [0.0 migration1to2; migration2to1 0.0]

    network = Network([population1, population2], migration_matrix)
    simulation = Simulation(network, dt, duration)

    try
        @time run_simulation!(simulation)
        thin_simulation!(simulation,thin_by) 
        
        open(joinpath(OUTPUT_DIRECTORY, "simulation_results_migration_$(migration_rate)_replicate_$(simulation_number).jld2"), "w") do file
            serialize(file, simulation)
        end
    catch e
        println("Error in simulation with args $args: $e")
    end
end

function main()
    # migration_rates = vcat([0.0], exp10.(LinRange(-6, 0.5, 9))) # Example migration rates to sweep over
    migration_rates = exp10.(LinRange(-6, 0.5, 9)) # Example migration rates to sweep over
    start_rep = 5001
    num_replicates = 10000

    # Creating a list of tuples with migration rates and simulation numbers
    simulation_args = [(rate, num) for rate in migration_rates for num in start_rep:(start_rep + num_replicates - 1)]

    @threads for args in simulation_args
        run_single_simulation(args)
    end
end

main()

