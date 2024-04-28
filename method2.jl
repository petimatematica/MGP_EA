## Projected Gradient Method with GPA1 ##

function GPA2(x, f, ∇f, projection, σ, min_step, β_inicial)
   
    β = β_inicial
    j = 0
 
    while β < min_step 
     z_kj = projection(x - β * 2.0^(-j) * ∇f(x)) # Calcula o ponto z_kj
     stptest = f(z_kj) - f(x) + σ * dot(∇f(x), x - z_kj) # Armijo
 
     #println(j, " ", stptest)
 
       if stptest > 0.0 # Armijo test   
        j += 1 
        else
        β = β_inicial * 2.0^(-j)
        return β
       end
    
    end
 
    println("Step length too small!")
 end

 function method2(f, ∇f, x, ε, max_iter, GPA2)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    iteration_time = Float64[]
    
    # Initialization
    iter = 0
    x = x0
    t0 = time()
    
    while true
        it0 = time()
        ∇fx = ∇f(x) 
        fx = f(x)
        gradnorm = norm(∇fx)
        
        # First stopping condition
        proj_stop = x - projection(x - ∇fx) 
        if norm(proj_stop) < ε
            println("The solution has found!")
            break
        end

        push!(fvals, fx)
        push!(gradnorms, gradnorm)
        
        # Update iteration
        iter += 1
        it = time() - it0   
        β = GPA2(x, f, ∇f, projection, σ, imax_iter, β_inicial)  
        x = projection(x - β * ∇fx) 
        push!(iteration_time, it)
        push!(stepsizes_γ, 1.0)
        push!(stepsizes_β, β) 
        
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