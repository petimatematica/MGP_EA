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
feasible_set = 2

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
   resultado = method1(f, ∇f, ε, max_iter, GPA1)
   elseif strategy == "GPA2"
   resultado = method2(f, ∇f, ε, max_iter, GPA2)
   elseif strategy == "GPA3"
   resultado = method3(f, ∇f, ε, max_iter) 
end

## Show the result ##

ENV["LINES"] = 1000
println(resultado[3])
# println("Real global minimum: ", global_minimum)
# println("Global minimum found: ", resultado[1])
println("Minimum value of f: ", resultado[2])
println("Total time spent: ", resultado[4])

## ORGANIZANDO OS TESTES ##

# projections = [projection1, projection2]
# dimensions = [2, 4]
# strategies = ["GPA1", "GPA2"]
#total_tests = length(feasible_sets) * length(dimensions) * length(strategies)
#random_x0 = [rand(dimensions[i]) for i in 1:2 for _ in 1:3 ]
#random_x0 = [[rand(dimension) for _ in 1:3] for dimension in dimensions]

# for k in 1:3
#    x0 = rand(4)
#    ENV["LINES"] = 10000
#    println(resultado[3])
#    println("Minimum value of f: ", resultado[2])
#    println("Total time spent: ", resultado[4])
#    println("x_0 = ", x0)   
# end 


















# for feasible_set in feasible_sets
#    if feasible_set == "projection1"
#       projection(x) = projection1(x)
#       elseif feasible_set == "projection2"
#       projection(x) = projection2(x)
#    end 
#    for strategy in strategies
#          if strategy == "GPA1"
#             resultado = method1(f, ∇f, ε, max_iter, GPA1)
#             elseif strategy == "GPA2"
#             resultado = method2(f, ∇f, ε, max_iter, GPA2)
#             elseif strategy == "GPA3"
#             resultado = method3(x0, f, ∇f, ε, max_iter)
#          end 
#       for dimension in dimensions
#          for k in 1:3
#              x0 = rand(dimension)
#              println("Executando teste com: conjunto viável = $feasible_set, dimensão = $dimension, estratégia = $strategy")
#              ENV["LINES"] = 10000
#              println(resultado[3])
#              println("Minimum value of f: ", resultado[2])
#              println("Total time spent: ", resultado[4])
#              println("x_0 = ", x0)   
#          end 
#       end 
#    end
# end

# Definir todas as estratégias, dimensões, projeções e chutes iniciais
# strategies = ["GPA1", "GPA2", "GPA3"]
# dimensions = [2, 3]  # Exemplo de duas dimensões
# feasible_sets = [1, 2]  # Exemplo de duas funções de projeção
# initial_guesses = [rand(d) for d in dimensions]  # Chutes iniciais aleatórios para cada dimensão

# for (x0_list, dimension) in zip(random_x0, dimensions)
#    for x0 in x0_list
#       for dimension in dimensions
#          for feasible_set in feasible_sets

#             if feasible_set == "1"
#             projection(x) = projection1(x)
#             elseif feasible_set == "2"
#             projection(x) = projection2(x)
#             elseif feasible_set == "3"
#             projection(x) = projection3(x)
#             elseif feasible_set == "4"
#             projection(x) = projection4(x)
#             elseif feasible_set == "5"
#             projection(x) = projection5(x)
#             elseif feasible_set == "6"
#             projection(x) = projection6(x)
#             end

#             for strategy in strategies
#                local resultado
#                if strategy == "GPA1"
#                resultado = method1(f, ∇f, ε, max_iter, GPA1)
#                elseif strategy == "GPA2"
#                resultado = method2(f, ∇f, ε, max_iter, GPA2)
#                elseif strategy == "GPA3"
#                resultado = method3(f, ∇f, ε, max_iter) 
#                end

#                println("Executando teste com: conjunto viável = $feasible_set, dimensão = $dimension, estratégia = $strategy")
#                ENV["LINES"] = 10000
#                println(resultado[3])
#                println("Minimum value of f: ", resultado[2])
#                println("Total time spent: ", resultado[4])
#                println("x_0", x0)
#             end
#          end
#       end
#    end
# end