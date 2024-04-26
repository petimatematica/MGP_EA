## Projected Gradient Method with GPA1 ##

function GPA2(x, f, ∇f, projection, σ, min_step, β_inicial)
   
    β = β_inicial
    j = 0
 
    while β < min_step 
     z_kj = projection(x - β * 2.0^(-j) * ∇f(x)) # Calcula o ponto z_kj
     stptest = f(z_kj) - f(x) + σ * dot(∇f(x), x - z_kj) # Testa Armijo para o z_kj obtido
 
     #println(j, " ", stptest)
 
       if stptest > 0.0 # Se a condição de Armijo não for satisfeita, testa o próximo j   
        j += 1 
        else
        β = β_inicial * 2.0^(-j)
        return β
       end
    
    end
 
    println("Step length too small!")
 end

 function method2(f, ∇f, x, ε, max_iter, GPA2)

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
        β = GPA2(x, f, ∇f, projection, σ, imax_iter, β_inicial)  # GPA2
        x = projection(x - β * ∇fx) # GPA2
        push!(stepsizes_γ, 1.0)
        push!(stepsizes_β, β) # Armazena o comprimento de passo na lista stepsizes      
    end

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ)
    
    return (x, f(x), info)
end