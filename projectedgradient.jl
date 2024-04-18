#
# Método do Gradiente Projetado
#
x = x0

function gradienteproj(f, ∇f, x, ε, max_iter, select_strategy)

    # Listas para armazenar os resultados
    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes = Float64[]
    
    # Inicialização
    # x = x0
    iter = 0
    # push!(stepsizes, NaN)
    # push!(fvals, f(x))
    # push!(gradnorms, ∇f(x))
    
    while true

        ∇fx = ∇f(x) # Calcula o gradiente no ponto atual
        fx = f(x) # Calcula o valor da função no ponto atual
        gradnorm = norm(∇fx) # Calcula a norma do gradiente no ponto atual
        
        # Verifica o critério de parada
        proj_stop = x - projection(x - ∇fx) 
        if norm(proj_stop) < ε || iter > max_iter
            break
        end

        push!(fvals, fx) # Armazena o valor da função na lista fvals
        push!(gradnorms, gradnorm) # Armazena o valor do gradiente da função na lista gradnorms
        
        # Atualiza o contador de iterações
        iter += 1

        # GPA1
            newsteps = GPA1(x, f, ∇f, projection, σ, imax_iter, γ_inicial) # GPA1
            z = projection(x - newsteps[2]  * ∇f(x)) # GPA1
            x = x + newsteps[1] * (z - x) # GPA1
            push!(stepsizes, newsteps[1])
        
        # GPA2
            # β = GPA2(x, f, ∇f, projection, σ, imax_iter, β_inicial)  # GPA2
            # x = projection(x - β * ∇fx) # GPA2
            # push!(stepsizes, β) # Armazena o comprimento de passo na lista stepsizes
    end

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes = stepsizes)
    
    return (x, f(x), info)
end

