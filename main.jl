## Main Code ##

include("GPA1.jl")
include("GPA2.jl")

using LinearAlgebra, DataFrames, Random, JLD2

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
   a = fill(-1, n)
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

guess = MersenneTwister(12345)
x_rand = rand(guess, 500)
n = length(x_rand) 
σ = 1.e-4 
ε = 1.e-5 
β_start = 1.0
β1 = 0.01
β2 = 0.9
γ_start = 1.0
min_step = 1.e-5
max_iter = 50000
strategy = "GPA1"
projection = projection6

# Definitions #

x0 = projection(x_rand)

if strategy == "GPA1"
   (x, f(x), info, et, ierror, seqx, evalsf) = method1(x0, f, ∇f, ε, max_iter, GPA1, projection)
   else
   (x, f(x), info, et, ierror, seqx, evalsf) = method2(x0, f, ∇f, ε, max_iter, GPA2, projection) 
end

# Show the result #

ENV["LINES"] = 10000
println(info)
println("Minimum value of f: ", f(x))
println("Total time spent: ", et)
println("Function evaluations = ", evalsf)
println("Ierror = ", ierror)
