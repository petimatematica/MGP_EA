## MAIN ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")
include("method3.jl")

using LinearAlgebra, DataFrames

## Parameters choose ##

x0 = [1, 1.9]
n = length(x0) # Dimensão da Função Dixon-Price
σ = 1.e-4 # Parâmetro da Busca de Armijo
ε = 1.e-5 # Critério de parada do Método
β_inicial = 1.0 # Comprimento de passo inicial
β1 = 1.e-6
β2 = 1.0
γ_inicial = 1.0
min_step = 1.e-5
max_iter = 10000
strategy = "GPA3"
feasible_set = "1"

## Estrategy choose ##

if strategy == "GPA1"
   resultado = method1(f, ∇f, ε, max_iter, GPA1)
   elseif strategy == "GPA2"
   resultado = method2(f, ∇f, ε, max_iter, GPA2)
   elseif strategy == "GPA3"
   resultado = method3(f, ∇f, ε, max_iter) 
end

## Feasible set choose ##

if feasible_set == "1"
   projection(x) = projection1(x)
   elseif projection == "2"
   projection(x) = projection2(x)
   elseif projection == "3"
   projection(x) = projection3(x)
   elseif projection == "4"
   projection(x) = projection4(x)
   elseif projection == "5"
   projection(x) = projection5(x)
   elseif projection == "6"
   projection(x) = projection6(x)
end

## Show the result ##

ENV["LINES"] = 1000
println(resultado[3])
println("Real global minimum: ", global_minimum)
println("Global minimum found: ", resultado[1])
println("Minimum value of f: ", resultado[2])
println("Total time spent: ", resultado[4])

## ORGANIZANDO OS TESTES ##

projection_functions = [projection1, projection2, projection3, projection4, projection5, projection6]
strategies = ["GPA1", "GPA2", "GPA3"]
ndims = [5, 10, 50, 100, 500]
initial_points = 