#
# Método do Gradiente Projetado
#

using LinearAlgebra, DataFrames

include("dixonpricedim2.jl")
include("projections.jl")
include("GPA2.jl")

function gradienteproj(f, ∇f, x0, max_step, η)
    
    # Parâmetros de controle
    ϵ = 1e-6 # Critério de parada para a norma do gradiente
    max_iter = 10000 # Número máximo de iteradas
    σ = 0.5


    # Listas para armazenar os resultados
    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes = Float64[]
    
    # Inicialização
    x = x0
    iter = 0
    
    while true

        ∇fx = ∇f(x) # Calcula o gradiente no ponto atual
        fx = f(x) # Calcula o valor da função no ponto atual
        gradnorm = norm(∇fx) # Calcula a norma do gradiente no ponto atual
        
        # Verifica o critério de parada
        proj_stop = x - projection(x - ∇fx) 
        if norm(proj_stop) < ϵ || iter >= max_iter
            break
        end

        push!(fvals, fx) # Armazena o valor da função na lista fvals
        push!(gradnorms, gradnorm) # Armazena o valor do gradiente da função na lista gradnorms
        
        # Busca de Armijo GPA2
        # β = max_step
        # while f(x - β * ∇fx) > f(x) - η * β * dot(∇fx, ∇fx) || α < minstep
        #     β *= 0.5
        # end

        push!(stepsizes, α) # Armazena o comprimento de passo na lista stepsizes
        
        # Atualiza o ponto
        # z = projection(x - β * ∇fx) # GPA1
        # x += γ * (z - x) # GPA1

        x = projection(x - β * ∇fx) # GPA2
        
        # Atualiza o contador de iterações
        iter += 1
    end

    info = DataFrame(fvals=fvals, gradnorms=gradnorms, stepsizes=stepsizes)
    
    return x, f(x), info
end

# Chamada da função gradiente_armijo com a função teste
x0 = [5; 1]  # Ponto inicial
max_step = 1.0    # Tamanho máximo do passo
η = 0.5           # Parâmetro de Armijo
minstep = 1e-4
resultado = gradiente_armijo(teste_funcao, grad_teste_funcao, x0, max_step, η)

# Exibir o resultado
println(resultado[3])
println("Ponto de mínimo: ", resultado[1])
println("Valor mínimo da função: ", resultado[2])

