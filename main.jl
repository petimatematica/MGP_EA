# Chamar todos os arquivos

include("dixonprice.jl")
include("projections.jl")
include("method1.jl")
include("method2.jl")
include("method3.jl")

using LinearAlgebra, DataFrames

# Definição dos Parâmetros

x0 = [1, 1.9, 2] # Ponto inicial
n = length(x0) # Dimensão da Função Dixon-Price
y = [1, 1,0] # Centro do conjunto viável (Bola)
δ = 1.0 # Raio da Bola
σ = 1.e-4 # Parâmetro da Busca de Armijo
ε = 1.e-4 # Critério de parada do Método
β_inicial = 1.0 # Comprimento de passo inicial
β1 = 1.e-6
β2 = 1.0
γ_inicial = 1.0
min_step = 1.e-6
max_iter = 100000 # Máximo de iteradas do Método
strategy = "GPA3"

if strategy == "GPA1"
   resultado = method1(f, ∇f, x0, ε, max_iter, GPA1) 
elseif strategy == "GPA2"
   resultado = method2(f, ∇f, x0, ε, max_iter, GPA2)
elseif strategy == "GPA3"
   resultado = method3(f, ∇f, x0, ε, max_iter) 
end

# Exibir o resultado
ENV["LINES"] = 1000
println(resultado[3])
println("Real global minimum: ", global_minimum)
println("Global minimum found: ", resultado[1])
println("Minimum value of f: ", resultado[2])
println("Time spent: ", resultado[4])

