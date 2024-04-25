# Chamar todos os arquivos

include("dixonprice.jl")
include("projections.jl")
include("GPA1.jl")
include("GPA2.jl")
include("projectedgradient.jl")

using LinearAlgebra, DataFrames

# Definição dos Parâmetros

x0 = [1, 1.9] # Ponto inicial
n = length(x0) # Dimensão da Função Dixon-Price
y = [1, 1] # Centro do conjunto viável (Bola)
δ = 1.0 # Raio da Bola
σ = 1.e-4 # Parâmetro da Busca de Armijo
ε = 1.e-5 # Critério de parada do Método
β_inicial = 1.0 # Comprimento de passo inicial
β1 = 1.e-6
β2 = 1.0
γ_inicial = 1.0
min_step = 1.e-6
max_iter = 500 # Máximo de iteradas do Método
imax_iter = 50 # Máximo de iteradas da estratégia

strategy = "GPA1"
resultado = gradienteproj(f, ∇f, x0, ε, max_iter, strategy)

# Exibir o resultado
ENV["LINES"] = 1000
println(resultado[3])
println("Ponto de mínimo: ", resultado[1])
println("Valor mínimo da função: ", resultado[2])
