#
# Gradient Project Method
#

using LinearAlgebra, DataFrames

function gradiente_armijo(f, ∇f, x0, max_step, η)
    
    # Parâmetros de controle
    ϵ = 1e-6 # Critério de parada para a norma do gradiente
    max_iter = 1000 # Número máximo de iteradas
    η =0.5

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
        if norm(∇fx) < ϵ || iter >= max_iter
            break
        end

        push!(fvals, fx) # Armazena o valor da função na lista fvals
        push!(gradnorms, gradnorm) # Armazena o valor do gradiente da função na lista gradnorms
        
        # Busca de Armijo
        α = max_step
        while f(x - α * ∇fx) > f(x) - η * α * dot(∇fx, ∇fx) || α < minstep
            α *= 0.5
        end

        push!(stepsizes, α) # Armazena o comprimento de passo na lista stepsizes
        
        # Atualiza o ponto
        x -= α * ∇fx
        
        # Atualiza o contador de iterações
        iter += 1
    end

    info = DataFrame(fvals=fvals, gradnorms=gradnorms, stepsizes=stepsizes)
    
    return x, f(x), info
end

# # Definição da função teste
function teste_funcao(x)
    return (x[1] - 1)^2 + 2 * (2 * x[2]^2 - x[1])^2
end

# Definição do gradiente da função teste
function grad_teste_funcao(x)
    return [6*x[1] - 8*x[2]^2 - 2, 32*x[2]^3 - 16*x[1]*x[2]]
end

# Chamada da função gradiente_armijo com a função teste

x0 = [1; 3]  # Ponto inicial
max_step = 1.0    # Tamanho máximo do passo
η = 0.5           # Parâmetro de Armijo
minstep = 1e-4
resultado = gradiente_armijo(teste_funcao, grad_teste_funcao, x0, max_step, η)

# Exibir o resultado
println(resultado[3])
println("Ponto de mínimo: ", resultado[1])
println("Valor mínimo da função: ", resultado[2])

