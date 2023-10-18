using Base.Threads
using Serialization
using Random
include("../coevolution_network_base.jl")
using .CoevolutionNetworkBase

const OUTPUT_DIRECTORY = "/pool001/dswartz/viral_coev_output/betasweep_results"

if !isdir(OUTPUT_DIRECTORY)
    mkdir(OUTPUT_DIRECTORY)
end

println("Number of threads: ", nthreads())

const L = 60.0
const dx = 0.3
const x = -L/2:dx:L/2-dx
const r = 3.5
const M = 15
const alpha = 1.0
const gamma = 0.0
const D = 0.02
const Nh = 3 * 10^6
const dt = 0.05
const duration = 80.0

# Find the index of the grid point closest to x = 0
index_closest_to_zero = argmin(abs.(x))

# Initialize viral_density as zeros
const viral_density = zeros(Float64, length(x))

# Set the value of viral_density at the closest index to 1/dx
viral_density[index_closest_to_zero] = 1/dx

const immune_density = zeros(Float64, length(x))

function run_single_simulation(args; thin_by=50)
    beta_value, simulation_number = args

    # Generating a unique seed using beta_value and simulation_number
    seed = hash((beta_value, simulation_number))
    Random.seed!(seed)

    population = Population(L, dx, r, M, beta_value, alpha, gamma, D, Nh, viral_density, immune_density)
    
    migration_matrix = [0.0][:,:]

    network = Network([population], migration_matrix)
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

job_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

function main(job_id)
    beta_values = LinRange(1.05, 6, 10)
    
    # Assign a specific beta value based on job_id
    beta_value = beta_values[job_id]
    
    num_replicates = 8000  # Total number of replicates for the assigned beta_value

    # Create a list of simulation numbers for the given beta_value
    simulation_numbers = 1:num_replicates

    @threads for simulation_number in simulation_numbers
        args = beta_value, simulation_number
        run_single_simulation(args)
    end
end

main(job_id)