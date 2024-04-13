# Chamar todos os arquivos

include("dixonpricedim2.jl")
include("projections.jl")
include("GPA2.jl")
include("projectedgradient.jl")

using LinearAlgebra, DataFrames

# Definição dos Parâmetros

n = 2 # Dimensão da Função Dixon-Price
y = [2,2] # Centro do conjunto viável (Bola)
δ = 4.0 # Raio da Bola
x0 = [1,2] # Ponto inicial
σ = 1.e-4 # Parâmetro da Busca de Armijo
ε = 1.e-5 # Critério de parada do Método
β_inicial = 1.0 # Comprimento de passo inicial
max_iter = 10000 # Máximo de iteradas do Método
imax_iter = 100 # Máximo de iteradas da estratégia

# Chamada da função gradientproj
resultado = gradienteproj(f, ∇f, x0, ε, max_iter, GPA2)

# Exibir o resultado
println(resultado[3])
println("Ponto de mínimo: ", resultado[1])
println("Valor mínimo da função: ", resultado[2])
