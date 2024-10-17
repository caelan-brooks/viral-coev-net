using Distributed
using Base.Threads
using Serialization
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const OUTPUT_DIRECTORY = "C:/Users/Daniel/Desktop/simresults_initial_variance"

if !isdir(OUTPUT_DIRECTORY)
    mkdir(OUTPUT_DIRECTORY)
end

println("Number of threads: ", nthreads())

const L = 40.0
const dx = 0.1
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const beta = 2.5
const alpha = 0.0
const gamma = 1.0
const D = 0.01
const Nh = 2 * 10^6
const dt = 0.05
const duration = 100.0
const sigma = 2.0
const noise_method = :PL_with_dx

# Initialize viral_density as zeros
const immune_density = zeros(Float64, length(x))

const thin_by = 10

using Random

function run_single_simulation(args)
    initial_variance, simulation_number = args

    # Generating a unique seed using beta_value and simulation_number
    seed = hash((initial_variance, simulation_number))
    Random.seed!(seed)

    # viral_density_unnorm = zeros(Float64, length(x))
    # viral_density_unnorm[abs.(x) .< initial_variance/2] .= 1
    # viral_density = 100 .* viral_density_unnorm ./ sum(viral_density_unnorm .* dx)
    xs = -L/2:dx:L/2-dx
    viral_density = 100 .* exp.(-xs.^2 ./ (2 * initial_variance)) ./ sqrt.(2 * pi * initial_variance)

    population1 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density; noise_method=noise_method, sigma=sigma)
    
    migration_matrix = [0.0][:,:]

    network = Network([population1], migration_matrix)
    simulation = Simulation(network, dt, duration; thin_by=thin_by);

    # @time run_simulation!(simulation)
    # thin_simulation!(simulation,thin_by) 
    
    # open(joinpath(OUTPUT_DIRECTORY, "simulation_results_beta_$(beta_value)_replicate_$(simulation_number).jld2"), "w") do file
    #     serialize(file, simulation)
    # end

    try
        @time run_simulation!(simulation)
        # thin_simulation!(simulation) 
        
        open(joinpath(OUTPUT_DIRECTORY, "simulation_results_variance_$(initial_variance)_replicate_$(simulation_number).jld2"), "w") do file
            serialize(file, simulation)
        end
    catch e
        println("Error in simulation with args $args: ", e)
        # println(stacktrace(catch_backtrace()))
    end
end

function main()
    variance_values = LinRange(0.01, 0.1, 6)  # variance values to sweep over
    start_rep = 1
    num_replicates = 3000

    # Creating a list of tuples with beta values and simulation numbers
    simulation_args = [(variance, num) for variance in variance_values for num in start_rep:(start_rep + num_replicates - 1)]

    @threads for args in simulation_args
        run_single_simulation(args)
    end
end

main()
