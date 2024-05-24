## MAIN ##

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")
include("method3.jl")

using LinearAlgebra, DataFrames, BenchmarkProfiles, Plots

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
strategies = ["GPA1", "GPA2", "GPA3"]
dimensions = [3, 8, 10]
nguess = 3

times = Float64[] # Lista para armazenar os tempos de execução
iters = Float64[]  # Lista para armazenar a quantidade de iteradas
avalf = Float64[] # Lista para armazenar a quantiddade de avaliações de função

for strategy in strategies
   for feasible_set in feasible_sets 
      for dimension in dimensions
         for k in 1:nguess
         global n
         x0 = rand(dimension)
         n = length(x0)
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
            elseif strategy == "GPA2"
            resultado = method2(x0, f, ∇f, ε, max_iter, GPA2)
            elseif strategy == "GPA3"
            resultado = method3(x0, f, ∇f, ε, max_iter, min_step) 
         end
         t1 = time()
         elapsed_time = t1 - t0
         
         if resultado[5] > 0
            push!(times, Inf)
            push!(iters, Inf)
            push!(avalf, Inf)
         else
            iteration = size(resultado[6], 2)
            total_avals = sum(resultado[7])
            push!(times, elapsed_time)
            push!(iters, iteration)
            push!(avalf, total_avals)
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

total = length(feasible_sets) * length(dimensions) * nguess
println("Total tests: ", 3*total)

## Performance profile ##

h = total;

X = [times[1:h] times[h+1:2h] times[2h+1:3h]]; #Matriz com os tempos
Y = [iters[1:h] iters[h+1:2h] times[2h+1:3h]]; #Matriz com as iteradas
Z = [avalf[1:h] avalf[h+1:2h] avalf[2h+1:3h]]; #Matriz com as avaliações de função

colors=[:royalblue1, :green2, :orange]

# P1 = performance_profile(PlotsBackend(), X, ["GPA1", "GPA2", "GPA3"], 
# xlabel = "CPU time ratio", ylabel = "Solved problems [%]", legend = :bottomright, 
# palette = colors, linewidth = 2)

P2 = performance_profile(PlotsBackend(), Y, ["GPA1", "GPA2", "GPA3"], 
xlabel = "Iteration", ylabel = "Solved problems [%]", legend = :bottomright, 
palette = colors, linewidth = 2.5)

# P3 = performance_profile(PlotsBackend(), Z, ["GPA1", "GPA2"], 
# xlabel = "Function evaluations", ylabel = "Solved problems [%]", legend = :bottomright, 
# palette = colors, linewidth = 2)

plot(P2)

# final = plot(P1, P2, P3, layout=(1,3), size=(1400, 300))
# savefig("Performance.png")

