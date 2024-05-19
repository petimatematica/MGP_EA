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
strategy = "GPA2"
feasible_set = 1

## Feasible set choose ##

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

## Estrategy choose ##

if strategy == "GPA1"
   resultado = method1(x0, f, ∇f, ε, max_iter, GPA1)
   elseif strategy == "GPA2"
   resultado = method2(x0, f, ∇f, ε, max_iter, GPA2)
   elseif strategy == "GPA3"
   resultado = method3(x0, f, ∇f, ε, max_iter) 
end

## Show the result ##

ENV["LINES"] = 1000
println(resultado[3])
# println("Global minimum found: ", resultado[1])
println("Minimum value of f: ", resultado[2])
println("Total time spent: ", resultado[4])
println("Encontrou a solução?", " ", resultado[5])