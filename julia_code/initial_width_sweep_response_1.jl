using Distributed
using Base.Threads
using Serialization
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase
using DataFrames
using CSV
using Statistics

const OUTPUT_DIRECTORY = "../final_plots_viral_coev/variances_and_probabilities_response.csv"

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

    try
        @time run_simulation!(simulation)
    catch e
        println("Error in simulation with args $args: ", e)
        # println(stacktrace(catch_backtrace()))
    end

    total_infected = calculate_total_infected(simulation)
    antigenic_variance = calculate_antigenic_variance_per_deme(simulation)

    max_infected_index = argmax(total_infected)
    variance_at_max = antigenic_variance[1,max_infected_index]

    if total_infected[end] > 0
        surv = 1
    else
        surv = 0
    end

    return variance_at_max, surv
end

function main()
    variance_values = LinRange(0.01, 0.1, 6)  # variance values to sweep over
    start_rep = 1
    num_replicates = 3000

    # Creating a list of tuples with beta values and simulation numbers
    simulation_args = [(variance, num) for variance in variance_values for num in start_rep:(start_rep + num_replicates - 1)]

    df = DataFrame(variance=Float64[], probability=Float64[], variance_error=Float64[], probability_error=Float64[])

    for (var_idx, var) in enumerate(variance_values)
        println("variance index: $(var_idx)")
        variance_replicates = fill(NaN, num_replicates)
        survival_replicates = fill(NaN, num_replicates)
        @threads for rep in start_rep:(start_rep + num_replicates - 1)
            var_at_max, surv = run_single_simulation((var, rep))

            variance_replicates[rep] = var_at_max
            survival_replicates[rep] = surv
        end
        avg_var_at_max = mean(variance_replicates)
        err_var_at_max = std(variance_replicates) / sqrt(num_replicates)
        avg_surv = mean(survival_replicates)
        err_surv = std(survival_replicates) / sqrt(num_replicates)
        push!(df, (avg_var_at_max, avg_surv, err_var_at_max, err_surv))
    end
    CSV.write("../final_plots_viral_coev/variances_and_probabilities_2.csv", df)
end

main()
