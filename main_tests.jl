## MAIN ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")
include("method3.jl")

using LinearAlgebra, DataFrames, ProgressMeter

## Parameters choose ##

σ = 1.e-4 # Parâmetro da Busca de Armijo
ε = 1.e-5 # Critério de parada do Método
β_inicial = 1.0 # Comprimento de passo inicial
β1 = 1.e-6
β2 = 1.0
γ_inicial = 1.0
min_step = 1.e-5
max_iter = 10000

## Organizando os testes ##

feasible_sets = [1]
strategies = ["GPA1", "GPA2"]
dimensions = [3, 10, 15]

total_iter = length(feasible_sets) * length(strategies) * 3
@showprogress 1 "Executando testes..." for iter in 1:total_iter

for dimension in dimensions 
   for feasible_set in feasible_sets 
      for strategy in strategies
         for k in 1:5
         global n
         x0 = rand(dimension)
         n = length(x0)
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

         if strategy == "GPA1"
            resultado = method1(x0, f, ∇f, ε, max_iter, GPA1)
            elseif strategy == "GPA2"
            resultado = method2(x0, f, ∇f, ε, max_iter, GPA2)
            elseif strategy == "GPA3"
            resultado = method3(x0, f, ∇f, ε, max_iter) 
         end

         ENV["LINES"] = 10000
         println(resultado[3])
         println("Minimum value of f: ", resultado[2])
         println("Total time spent: ", resultado[4])
         println("x_0 = ", x0) 
         println("Executando teste com: conjunto viável = $feasible_set, estratégia = $strategy, dimension = $dimension")  
         end
      end
   end 
end
end

total = length(strategies) * 3 * length(feasible_sets) * length(dimensions)
println("Total problems: ", total)