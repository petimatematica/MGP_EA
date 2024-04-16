# Chamar todos os arquivos

include("dixonprice.jl")
include("projections.jl")
include("GPA2.jl")
include("projectedgradient.jl")

using LinearAlgebra, DataFrames

# Definição dos Parâmetros

x = [1,2] # Ponto inicial
n = length(x) # Dimensão da Função Dixon-Price
y = [1,0] # Centro do conjunto viável (Bola)
δ = 5.0 # Raio da Bola
σ = 1.e-4 # Parâmetro da Busca de Armijo
ε = 1.e-5 # Critério de parada do Método
β_inicial = 1.0 # Comprimento de passo inicial
max_iter = 10000 # Máximo de iteradas do Método
imax_iter = 100 # Máximo de iteradas da estratégia

# Chamada da função gradientproj
resultado = gradienteproj(f, ∇f, x, ε, max_iter, GPA2)

# Exibir o resultado
ENV["LINES"] = 1000
println(resultado[3])
println("Ponto de mínimo: ", resultado[1])
println("Valor mínimo da função: ", resultado[2])
