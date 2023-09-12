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
ENV["JULIA_NUM_THREADS"] = 12

println("Number of threads: ", nthreads())

# Define function to run a single simulation
function run_single_simulation(args)
    println(args)
    migration_rate, simulation_number = args

    # The parameters and initializations here are the same as in your original script
    L = 40.0
    dx = 0.3
    x = -L/2:dx:L/2-dx
    r = 3.0
    M = 15
    beta = 2.0
    alpha = 1.0
    gamma = 0.0
    D = 0.01
    Nh = 1 * 10^6
    dt = 0.05
    duration = 80.0

    viral_density = [abs(val) <= 0.5 ? 100.0 : 0.0 for val in x]
    viral_density2 = zeros(Float64, length(x))
    immune_density = zeros(Float64, length(x))

    population1 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density)
    population2 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density2, immune_density)
    
    migration_matrix = [0 migration_rate; migration_rate 0]

    network = Network([population1, population2], migration_matrix)
    simulation = Simulation(network, dt, duration)

    @time begin
        run_simulation!(simulation)
    end

    total_infected = calculate_total_infected(simulation)

    result_dict = Dict("times" => simulation.duration_times, "total_infected_number" => total_infected)
    
    open("simresults_largedt/simulation_results_migration_$(migration_rate)_replicate_$(simulation_number).jld2", "w") do file
        serialize(file, result_dict)
    end
end

function main()
    migration_rates = exp10.(LinRange(-6, 0.5, 9)) # Example migration rates to sweep over
    start_rep = 0
    num_replicates = 60

    # Creating a list of tuples with migration rates and simulation numbers
    simulation_args = [(rate, num) for rate in migration_rates for num in start_rep:(start_rep + num_replicates - 1)]

    @threads for args in simulation_args
        run_single_simulation(args)
    end
end

main()
