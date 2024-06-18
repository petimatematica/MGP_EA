## Main Code ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")

using LinearAlgebra, DataFrames, Random

# Parameters #
guess = MersenneTwister(1234)
x = rand(guess, 50)
n = length(x) 
σ = 1.e-4 
ε = 1.e-5 
β_start = 1.0
β1 = 1.e-6
β2 = 1.0
γ_start = 1.0
min_step = 1.e-5
max_iter = 30000
strategy = "GPA2"
feasible_set = 6

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

x0 = projection(x)

if strategy == "GPA1"
   result = method1(x0, f, ∇f, ε, max_iter, GPA1)
   elseif strategy == "GPA2"
   result = method2(x0, f, ∇f, ε, max_iter, GPA2)
   elseif strategy == "GPA3"
   result = method3(x0, f, ∇f, ε, max_iter, min_step) 
end

# Show the result #

ENV["LINES"] = 1000
println(result[3])
println("Minimum value of f: ", result[2])
println("Total time spent: ", result[4])
println("Function evaluations = ", sum(result[7]))
println("Ierror = ", result[5])
