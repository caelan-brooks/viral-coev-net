using Plots
using LaTeXStrings
using Statistics
using Base.Threads
include("../coevolution_network_base.jl")
using .CoevolutionNetworkBase
using Printf

# Parameters
L = 40.0
dx = 0.05
x = -L/2:dx:L/2-dx
r = 3.0
M = 15
beta = 2.5
alpha = 0.0
gamma = 1.0
# D = 0.0025
D = 0.01
Nh = 2 * 10^6
stochastic = true
sigma = 1.0
dt = 0.05
duration = 20.0

migration_rates = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0]
n_replicates = 200
colors = distinguishable_colors(length(migration_rates))

FSTs_mean = Vector{Vector{Float64}}()
FSTs_std = Vector{Vector{Float64}}()
ts_ref = nothing

function run_replicate(mig_rate)
    viral_density = zeros(length(x))
    viral_density[Int(round(length(x)/2))] = 100/dx
    viral_density2 = zeros(length(x))
    immune_density = zeros(length(x))

    pop1 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density, immune_density; stochastic=stochastic, sigma=sigma)
    pop2 = Population(L, dx, r, M, beta, alpha, gamma, D, Nh, viral_density2, immune_density; stochastic=stochastic, sigma=sigma)
    pops = [pop1, pop2]

    mig_matrix = mig_rate * ones(2, 2)
    net = Network(pops, mig_matrix)
    sim = Simulation(net, dt, duration; thin_by=1)
    run_simulation!(sim)

    return calculate_FST(sim), sim.duration_times
end

for (i, mig_rate) in enumerate(migration_rates)
    println("Running migration rate $mig_rate")
    FST_reps = Vector{Vector{Float64}}(undef, n_replicates)

    Threads.@threads for rep in 1:n_replicates
        FST_reps[rep], ts_local = run_replicate(mig_rate)
    end


    FST_mat = hcat(FST_reps...)  # [time × replicate]
    push!(FSTs_mean, mean(FST_mat, dims=2)[:])
    push!(FSTs_std, std(FST_mat, dims=2)[:])
end

_, ts_ref = run_replicate(0.0)

using CSV, DataFrames

for (i, mean_val) in enumerate(FSTs_mean)
    df = DataFrame(time = ts_ref, FST = mean_val, FST_error = FSTs_std[i])
    filename = "FST_m$(migration_rates[i]).csv"
    CSV.write(filename, df)
end
