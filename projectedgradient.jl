#
# Método do Gradiente Projetado
#
#x = x0

function gradienteproj(f, ∇f, x, ε, max_iter, strategy)

    # Listas para armazenar os resultados
    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    
    # Inicialização
    iter = 0
    x = x0
    
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

        if strategy == "GPA1"
        # GPA1
            newsteps = GPA1(x, f, ∇f, projection, σ, imax_iter, γ_inicial, β_inicial) # GPA1
            z = projection(x - newsteps[2]  * ∇f(x)) # GPA1
            x = x + newsteps[1] * (z - x) # GPA1
            push!(stepsizes_β, newsteps[2])
            push!(stepsizes_γ, newsteps[1])   
        else
        # GPA2
            β = GPA2(x, f, ∇f, projection, σ, imax_iter, β_inicial)  # GPA2
            x = projection(x - β * ∇fx) # GPA2
            push!(stepsizes_γ, 1.0)
            push!(stepsizes_β, β) # Armazena o comprimento de passo na lista stepsizes     
        end    
    end

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ)
    
    return (x, f(x), info)
end

