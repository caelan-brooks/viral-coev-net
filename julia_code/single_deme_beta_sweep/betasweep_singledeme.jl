using Distributed
using Base.Threads
using Serialization
include("../coevolution_network_base.jl")
using .CoevolutionNetworkBase

const OUTPUT_DIRECTORY = "C:/Users/Daniel/Desktop/simresults_single_deme_beta_sweep"

if !isdir(OUTPUT_DIRECTORY)
    mkdir(OUTPUT_DIRECTORY)
end

println("Number of threads: ", nthreads())

const L = 50.0
const dx = 0.3
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const alpha = 1.0
const gamma = 0.0
const D = 0.01
const Nh = 3 * 10^6
const dt = 0.05
const duration = 70.0

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
    beta_value, simulation_number = args

    # Generating a unique seed using beta_value and simulation_number
    seed = hash((beta_value, simulation_number))
    Random.seed!(seed)

    population1 = Population(L, dx, r, M, beta_value, alpha, gamma, D, Nh, viral_density, immune_density)
    
    migration_matrix = [0.0][:,:] # this horrible line makes it a 1x1 matrix...

    network = Network([population1], migration_matrix)
    simulation = Simulation(network, dt, duration)

    try
        @time run_simulation!(simulation)
        thin_simulation!(simulation,thin_by) 
        
        open(joinpath(OUTPUT_DIRECTORY, "simulation_results_beta_$(beta_value)_replicate_$(simulation_number).jld2"), "w") do file
            serialize(file, simulation)
        end
    catch e
        println("Error in simulation with args $args: ", e)
        println(stacktrace(catch_backtrace()))
    end
end

function main()
    beta_values = LinRange(1.1, 6, 10)  # Example beta values to sweep over
    start_rep = 0
    num_replicates = 10000

    # Creating a list of tuples with beta values and simulation numbers
    simulation_args = [(beta, num) for beta in beta_values for num in start_rep:(start_rep + num_replicates - 1)]

    @threads for args in simulation_args
        run_single_simulation(args)
    end
end

main()
