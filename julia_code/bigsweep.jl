using DelimitedFiles
using Distributed
using Base.Threads
using Serialization
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

if !isdir("simresults_largedt")
    mkdir("simresults_largedt")
end

using Base.Threads

println("Number of threads: ", nthreads())

const L = 40.0
const dx = 0.3
const x = -L/2:dx:L/2-dx
const r = 3.0
const M = 15
const beta = 2.0
const alpha = 1.0
const gamma = 0.0
const D = 0.01
const Nh = 1 * 10^6
const dt = 0.05
const duration = 80.0

const viral_density = [abs(val) <= 0.5 ? 100.0 : 0.0 for val in x]
const viral_density2 = zeros(Float64, length(x))
const immune_density = zeros(Float64, length(x))

function run_single_simulation(args)
    migration_rate, simulation_number = args

    population1 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density)
    population2 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density2, immune_density)
    
    migration_matrix = [0 migration_rate; migration_rate 0]

    network = Network([population1, population2], migration_matrix)
    simulation = Simulation(network, dt, duration)

    @time run_simulation!(simulation)

    total_infected = calculate_total_infected(simulation)

    result_dict = Dict("times" => simulation.duration_times, "total_infected_number" => total_infected)
    
    open("simresults_largedt/simulation_results_migration_$(migration_rate)_replicate_$(simulation_number).jld2", "w") do file
        serialize(file, result_dict)
    end
end

function main()
    migration_rates = exp10.(LinRange(-6, 0.5, 9)) # Example migration rates to sweep over
    start_rep = 0
    num_replicates = 6000

    # Creating a list of tuples with migration rates and simulation numbers
    simulation_args = [(rate, num) for rate in migration_rates for num in start_rep:(start_rep + num_replicates - 1)]

    @threads for args in simulation_args
        run_single_simulation(args)
    end
end

main()
