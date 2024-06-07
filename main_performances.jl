## MAIN ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")
include("method3.jl")

using LinearAlgebra, DataFrames, BenchmarkProfiles, Plots, Random

## Parameters choose ##

σ = 1.e-4 # Parâmetro da Busca de Armijo
ε = 1.e-5 # Critério de parada do Método
β_start = 1.0 # Comprimento de passo inicial
β1 = 1.e-6
β2 = 1.0
γ_start = 1.0
min_step = 1.e-5
max_iter = 30000

## Organizando os testes ##
## Avaliando tempo, número de iteradas e número de avaliação de função ##

feasible_sets = [1]
strategies = ["GPA1", "GPA2"]
dimensions = [3, 5, 8]
guess = MersenneTwister(1234)
nguess = 3

times1 = Float64[] # Lista para armazenar os tempos de execução da GPA1
times2 = Float64[] # Lista para armazenar os tempos de execução da GPA2
iters1 = Float64[]  # Lista para armazenar a quantidade de iteradas da GPA1
iters2 = Float64[]  # Lista para armazenar a quantidade de iteradas da GPA2
avalf1 = Float64[] # Lista para armazenar a quantiddade de avaliações de função da GPA1
avalf2 = Float64[] # Lista para armazenar a quantiddade de avaliações de função da GPA2

for k in 1:nguess
   for feasible_set in feasible_sets 
      for dimension in dimensions
         global n
         x0 = rand(guess, dimension)
         n = length(x0)

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

            ENV["LINES"] = 10000
            println(resultado[3])
            println("Minimum value of f: ", resultado[2])
            println("Total time spent: ", resultado[4])
            println("x_0 = ", x0) 
            println("Executando teste com: conjunto viável = $feasible_set, estratégia = $strategy, dimension = $dimension")
            println("Iters = ", iters)
            println("Elapsed time = ", times)
            println("Function evaluations =", avalf)
         end
      end
   end 
end

total = length(feasible_sets) * length(dimensions) * nguess * 2
println("Total tests: ", total)

## Performance profile ##

X = [times1 times2]; #Matriz com os tempos
Y = [iters1 iters2]; #Matriz com as iteradas
Z = [avalf1 avalf2]; #Matriz com as avaliações de função

colors=[:blue2, :orangered2]

Plots.backend(:png)

P1 = performance_profile(PlotsBackend(), X, ["GPA1", "GPA2"], 
xlabel = "CPU time ratio", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2.5)

P2 = performance_profile(PlotsBackend(), Y, ["GPA1", "GPA2"], 
xlabel = "Iteration", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2.5)

P3 = performance_profile(PlotsBackend(), Z, ["GPA1", "GPA2"], 
xlabel = "Function evaluations", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2.5)

final = plot(P1, P2, P3, layout = (3, 1), size=(600,1500), left_margin = 10Plots.mm, dpi=300)

savefig(final, "performances_profile2.png")

