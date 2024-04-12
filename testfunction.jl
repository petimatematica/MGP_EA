#
#   Test function
#

using LinearAlgebra, DataFrames, Random, Printf

include("armijo.jl")
include("gradient.jl")

# function f(x)
#     return (x[1] - 1)^2 + 2 * (2 * x[2]^2 - x[1])^2
# end

# function ∇(x)
#     df_dx = 6*x[1] - 8*x[2]^2 - 2
#     df_dy = 32*x[2]^3 - 16*x[1]*x[2]
#     return [df_dx, df_dy]
# end

# f(x) = (x - 2)^2  # Função objetivo simples: (x - 2)²
# ∇f(x) = 2 * (x - 2)  # Gradiente da função objetivo

# Parâmetros de otimização
x0 = 5000  # Ponto inicial
ϵ = 1e-6  # Tolerância de convergência
η = 0.5   # Parâmetro de busca de Armijo
maxiter = 1000  # Número máximo de iterações
minstep = 1e-6  # Tamanho mínimo do passo

# Chamada da função de gradiente com busca de Armijo
(x_min, ierror, info, et, seqx) = gradiente(x0, f, ∇f, ϵ, η, maxiter, minstep, armijo)

# Exibindo os resultados
println("Resultado da otimização:")
println("Ponto de mínimo: $x_min")
println("Código de erro: $ierror")
println("Tempo decorrido: $et segundos")
println("Sequência de pontos visitados:")
show(seqx)

