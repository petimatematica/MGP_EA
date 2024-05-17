## Projected Gradient Method with GPA1 ##

function GPA2(x, f, ∇f, projection, σ, min_step, β_inicial)
   
    β = β_inicial
    j = 0
 
    while β > min_step 
     z_kj = projection(x - β * 2.0^(-j) * ∇f(x)) # Calcula o ponto z_kj
     stptest = f(z_kj) - f(x) + σ * dot(∇f(x), x - z_kj) # Armijo
 
     #println(j, " ", stptest)
 
       if stptest > 0.0 # Armijo test   
        j += 1 
        else
        β = β_inicial * 2.0^(-j)
        break
       end
    end
    
    return β
 
    println("Step length too small!")
 end

 function method2(x0, f, ∇f, ε, max_iter, GPA2)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    iteration_time = Float64[]
    
    # Initialization
    iter = 0
    x = x0
    t0 = time()

    if norm(x - projection(x - ∇f(x))) < ε
        println("x0 is a stationary point!")
        return (x, f(x))
    end
    
    while true
        x_k = copy(x)
        it0 = time()
        ∇fx = ∇f(x) 
        fx = f(x)
        gradnorm = norm(∇fx)
        β = GPA2(x, f, ∇f, projection, σ, min_step, β_inicial)  
        x = projection(x - β * ∇fx)
        it = time() - it0
        push!(fvals, fx)
        push!(gradnorms, gradnorm)
        push!(iteration_time, it)
        push!(stepsizes_γ, 1.0)
        push!(stepsizes_β, β) 
        
        # First stopping condition
        if norm(x - x_k) < ε
            println("The solution has found!")
            break
        end
        
        # Update iteration
        iter += 1

        # Second stopping condition
        if iter > max_iter
            println("Maximum of iterations was achieved! Stopping...")
            break
        end   
    end

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ, iteration_time = iteration_time)
    et = time() - t0
    return (x, f(x), info, et)
end