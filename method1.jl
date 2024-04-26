## Projected Gradient Method with GPA1 ##

function GPA1(x, f, ∇f, projection, σ, min_step, γ_inicial, β_inicial)
   
    β = β_inicial
    γ = γ_inicial
    j = 0
    z_k = projection(x - β * ∇f(x)) # Calcula o ponto z_k

    while γ < min_step 
     stptest = f(x + 2.0^(-j) * (z_k - x)) - f(x) - σ * 2.0^(-j) * dot(∇f(x), z_k - x) # Testa Armijo para o z_k obtido 
       #println(j, " ", stptest)
 
       if stptest > 0.0 # Se a condição de Armijo não for satisfeita, testa o próximo j   
        j += 1 
        else
        γ = 2.0^(-j)
        break
       end       
    end
    β = (- dot(z_k - x, ∇f(x)) * β^2) / (2.0 * (f(x + β * (z_k - x)) - f(x) - β * dot(z_k - x, ∇f(x))))
    println("β = ",β)
    if β < β1 || β > β2
       β = β / 2.0
       #println("β = ", " ", β)
    end
    return (γ, β)
 
    println("Step length too small!.")
 end

 function method1(f, ∇f, x, ε, max_iter, GPA1)

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
        newsteps = GPA1(x, f, ∇f, projection, σ, imax_iter, γ_inicial, β_inicial) # GPA1
        z = projection(x - newsteps[2]  * ∇f(x)) # GPA1
        x = x + newsteps[1] * (z - x) # GPA1
        push!(stepsizes_β, newsteps[2])
        push!(stepsizes_γ, newsteps[1])        
    end

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ)
    
    return (x, f(x), info)
end