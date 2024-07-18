## Main Code - Performance profiles ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")

using LinearAlgebra, DataFrames, BenchmarkProfiles, Plots, Random, JLD2

## Parameters ##

σ = 1.e-4 
ε = 1.e-5 
β_start = 1.0 
β1 = 1.e-6
β2 = 1.0
γ_start = 1.0
min_step = 1.e-5
max_iter = 50000

## Analysis of time, number of iterations and number of function evaluations ##

projections = [projection1, projection2]
strategies = ["GPA1", "GPA2"]
dimensions = [5, 10]
guess = MersenneTwister(1234)
nguess = 2

times1 = Float64[] 
times2 = Float64[] 
iters1 = Float64[]  
iters2 = Float64[]  
evalf1 = Float64[] 
evalf2 = Float64[]

t_start = time()

for k in 1:nguess
   for dimension in dimensions
      for projection in projections

         x_start = rand(guess, dimension)
         x0 = projection(x_start)
         global f
         global x0

         for strategy in strategies
         
            t0 = time()
            println("Running test with: projection = $projection, strategy = $strategy, dimension = $dimension")

            if strategy == "GPA1"
               (x, f(x), info, et, ierror, seqx, evalsf) = method1(x0, f, ∇f, ε, max_iter, GPA1, projection)
               t1 = time()
               elapsed_time = t1 - t0

               if ierror > 0
                  push!(times1, Inf)
                  push!(iters1, Inf)
                  push!(evalf1, Inf)
               else
                  iteration = size(seqx, 2)
                  total_evals = evalsf
                  push!(times1, elapsed_time)
                  push!(iters1, iteration)
                  push!(evalf1, total_evals)
               end 
               times = times1
               iters = iters1
               evalf = evalf1
               
            else
               (x, f(x), info, et, ierror, seqx, evalsf) = method2(x0, f, ∇f, ε, max_iter, GPA2, projection)
               t1 = time()
               elapsed_time = t1 - t0

               if ierror > 0
                  push!(times2, Inf)
                  push!(iters2, Inf)
                  push!(evalf2, Inf)
               else
                  iteration = size(seqx, 2)
                  total_evals = sum(evalsf)
                  push!(times2, elapsed_time)
                  push!(iters2, iteration)
                  push!(evalf2, total_evals)
               end 
               times = times2
               iters = iters2
               evalf = evalf2

            end
            filename = "echo/" * string("guess", k) * string("set", feasible_set) * string("dim", dimension) * strategy * string("ierror", ierror) * ".jld2"
            @save filename info
         end
      end
   end 
end

t_final = time() - t_start
println("Total time spent = ", t_final/60, " ", "minutes")
problems = length(projections) * length(dimensions) * nguess
println("Number of problems: ", problems)
println("Problems tested, generating performance profiles...")

## Performance profiles ##

X = [times1 times2]; 
Y = [iters1 iters2]; 
Z = [evalf1 evalf2]; 

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