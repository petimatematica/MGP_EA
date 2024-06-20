## Main Code - Performance profiles ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")

using LinearAlgebra, DataFrames, BenchmarkProfiles, Plots, Random

## Parameters choose ##

σ = 1.e-4 
ε = 1.e-5 
β_start = 1.0 
β1 = 1.e-6
β2 = 1.0
γ_start = 1.0
min_step = 1.e-5
max_iter = 30000

## Analysis of time, number of iterations and number of function evaluations ##

feasible_sets = [1, 2, 3, 4, 5, 6]
strategies = ["GPA1", "GPA2"]
dimensions = [5, 10, 50, 100, 500]
guess = MersenneTwister(1234)
nguess = 5

times1 = Float64[] 
times2 = Float64[] 
iters1 = Float64[]  
iters2 = Float64[]  
avalf1 = Float64[] 
avalf2 = Float64[]

t_inicial = time()

for k in 1:nguess
   for feasible_set in feasible_sets 
      for dimension in dimensions
         global n
         x = rand(guess, dimension)
         n = length(x)

         for strategy in strategies
         
            t0 = time()
            global projection

            if feasible_set == 1
               projection = projection1
               elseif feasible_set == 2
               projection = projection2
               elseif feasible_set == 3
               projection = projection3
               elseif feasible_set == 4
               projection = projection4
               elseif feasible_set == 5
               projection = projection5
               elseif feasible_set == 6
               projection = projection6
            end

            x0 = projection(x)
            println("Running test with: feasible set = $feasible_set, strategy = $strategy, dimension = $dimension")

            if strategy == "GPA1"
               resultado = method1(x0, f, ∇f, ε, max_iter, GPA1)
               t1 = time()
               elapsed_time = t1 - t0

               if resultado[5] > 0
                  push!(times1, Inf)
                  push!(iters1, Inf)
                  push!(avalf1, Inf)
               else
                  iteration = size(resultado[6], 2)
                  total_avals = sum(resultado[7])
                  push!(times1, elapsed_time)
                  push!(iters1, iteration)
                  push!(avalf1, total_avals)
               end 
               times = times1
               iters = iters1
               avalf = avalf1
               
            else
               resultado = method2(x0, f, ∇f, ε, max_iter, GPA2)
               t1 = time()
               elapsed_time = t1 - t0

               if resultado[5] > 0
                  push!(times2, Inf)
                  push!(iters2, Inf)
                  push!(avalf2, Inf)
               else
                  iteration = size(resultado[6], 2)
                  total_avals = sum(resultado[7])
                  push!(times2, elapsed_time)
                  push!(iters2, iteration)
                  push!(avalf2, total_avals)
               end 
               times = times2
               iters = iters2
               avalf = avalf2
            end
         end
      end
   end 
end

t_final = time() - t_inicial
println("Total time spent = ", t_final/60)
problems = length(feasible_sets) * length(dimensions) * nguess
println("Number of problems: ", problems)
println("Problems tested, generating performance profiles...")

## Performance profiles ##

X = [times1 times2]; 
Y = [iters1 iters2]; 
Z = [avalf1 avalf2]; 

colors=[:blue2, :green2]

P1 = performance_profile(PlotsBackend(), X, ["GPA1", "GPA2"], 
xlabel = "CPU time ratio", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2, dpi = 1000)

P2 = performance_profile(PlotsBackend(), Y, ["GPA1", "GPA2"], 
xlabel = "Iteration", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2, dpi = 1000)

P3 = performance_profile(PlotsBackend(), Z, ["GPA1", "GPA2"], 
xlabel = "Function evaluations", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2, dpi = 1000)

println("Performance profiles generated, saving figures...")

savefig(P1, "performance_profile_CPU_time_ratio.png")
savefig(P2, "performance_profile_iteration.png")
savefig(P3, "performance_profile_function_evaluations.png")

println("Figures saved, the code has finished running.")

