## Main Code ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")

using LinearAlgebra, DataFrames, Random

# Parameters #
guess = MersenneTwister(1234)
x_rand = rand(guess, 2)
n = length(x_rand) 
σ = 1.e-4 
ε = 1.e-5 
β_start = 1.0
β1 = 1.e-6
β2 = 1.0
γ_start = 1.0
min_step = 1.e-5
max_iter = 30000
strategy = "GPA1"
feasible_set = 1

# Conditions #

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

x0 = projection(x_rand)

if strategy == "GPA1"
   (x, f(x), info, et, ierror, seqx, evalsf) = method1(x0, f, ∇f, ε, max_iter, GPA1)
   elseif strategy == "GPA2"
   (x, f(x), info, et, ierror, seqx, evalsf) = method2(x0, f, ∇f, ε, max_iter, GPA2) 
end

# Show the result #

ENV["LINES"] = 1000
println(info)
println("Minimum value of f: ", f(x))
println("Total time spent: ", et)
println("Function evaluations = ", evalsf)
println("Ierror = ", ierror)
