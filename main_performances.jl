## Main Code - Performance profiles ##

include("GPA1.jl")
include("GPA2.jl")

using LinearAlgebra, DataFrames, BenchmarkProfiles, Plots, Random, JLD2

# Dixon-Price Function #

function f(x::Vector) 
   n = length(x) 
   sum = (x[1] - 1)^2  
   for i in 2:n  
       sum += i * (2*x[i]^2 - x[i-1])^2 
   end
   return sum 
end

function ∇f(x::Vector) 
   n = length(x)   
   g = zeros(n)    
   g[1] = 2 * (x[1] - 1) - 4 * (2 * x[2]^2 - x[1])  
   for i in 2:(n-1)
       g[i] = 2 * i * (2 * x[i]^2 - x[i-1]) * (4 * x[i]) - (i+1) * 2 * (2 * x[i+1]^2 - x[i])  
   end
   g[n] = n * 2 * (2 * x[n]^2-x[n-1]) * (4 * x[n])
   return g
end

#
# Projection functions for diferent feasible sets 
#

# Projection 01: Ball of ray δ and center y #

function projection1(x)
   n = length(x)
   c = fill(0, n) 
   δ = 0.5 
   cx = norm(c-x)
   if cx <= δ
       return x
   else
       return c - (δ*(c-x)/cx)
   end
end

# Projection 02: C is R^n+ #

function projection2(x)
   return max.(0, x)
end

# Projection 03: C is a hypercube of dimension n. For each dimension, the range is [-1,1] #

function projection3(x)
   n = length(x)
   a = fill(-1, n)
   b = fill(1, n)
   for j in 1:n   
       if x[j] < a[j]
          x[j] = a[j]
       elseif a[j] <= x[j] <= b[j]
          x[j] = x[j]
       elseif x[j] > b[j]
          x[j] = b[j] 
       end
   end
   return x
end

# Projection 04: C is the set of points that stisfy an equality with inner product #

function projection4(x)
   n = length(x)
   a = fill(1, n)
   b = 0
   return x + ((b - dot(a, x))/norm(a)^2)*a
end

# Projection 05: C is the set of points that stisfy an equality with inner product #

function projection5(x)
   n = length(x)
   a = fill(1, n)
   b = 0
   return x + ((min(0, b - dot(a, x)))/norm(a)^2)*a
end

# Projection 06: C is the set of points that satisfy an equality with a product between a vector and a matrix with defined rank #

function matrix(rank, n)
   rng = MersenneTwister(3214)
   matrix_0 = rand(rng, rank, n)  
   U, Σ, V = svd(matrix_0)
   Σ[rank + 1:end] .= 0      
   A = U[:, 1:rank] * Diagonal(Σ[1:rank]) * V[:, 1:rank]'
   return A
end

function projection6(x)
   n = length(x)
   rank = 2  
   A = matrix(rank, n)
   b = fill(0, size(A, 1))
   return x - A' * inv(A * A') * (A * x - b)
end

# Parameters #

σ = 1.e-4 
ε = 1.e-5 
β_start = 1.0 
β1 = 0.01
β2 = 0.9
γ_start = 1.0
min_step = 1.e-5
max_iter = 40000

## Analysis of time, number of iterations and number of function evaluations ##

projections = [projection1, projection2, projection3, projection5]
#, projection4, projection5, projection6]
strategies = ["GPA1", "GPA2"]
dimensions = [5, 10, 30, 50, 100, 150, 200, 300, 400, 500]
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

               # ENV["LINES"] = 10000
               # println(info)
               # println("Minimum value of f: ", f(x))
               # println("Total time spent: ", et)
               # println("Function evaluations = ", evalsf)
               # println("Ierror = ", ierror)

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

               ENV["LINES"] = 1000
               println(info)
               println("Minimum value of f: ", f(x))
               println("Total time spent: ", et)
               println("Function evaluations = ", evalsf)
               println("Ierror = ", ierror)

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

            filename = "echo/" * string("guess", k) * string(projection) * string("dim", dimension) * strategy * string("ierror", ierror) * string("min", f(x)) * ".jld2"
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

# Performance profiles ##

X = [times1 times2]; 
Y = [iters1 iters2]; 
Z = [evalf1 evalf2]; 

colors=[:blue2, :green2]

P1 = performance_profile(PlotsBackend(), X, ["GPA1", "GPA2"], 
xlabel = "CPU time ratio", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 1.5, dpi = 1000)

P2 = performance_profile(PlotsBackend(), Y, ["GPA1", "GPA2"], 
xlabel = "Iteration", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 1.5, dpi = 1000)

P3 = performance_profile(PlotsBackend(), Z, ["GPA1", "GPA2"], 
xlabel = "Function evaluations", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 1.5
, dpi = 1000)

println("Performance profiles generated, saving figures...")

savefig(P1, "performance_profile_CPU_time_ratio.png")
savefig(P2, "performance_profile_iteration.png")
savefig(P3, "performance_profile_function_evaluations.png")

println("Figures saved, the code has finished running.")