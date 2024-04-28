## Projected Gradient Method with GPA1 ##

function GPA1(x, f, ∇f, projection, σ, min_step, γ_inicial, β_inicial)
   
    β = β_inicial
    γ = γ_inicial
    j = 0
    z_k = projection(x - β * ∇f(x))

    while γ < min_step 
     stptest = f(x + 2.0^(-j) * (z_k - x)) - f(x) - σ * 2.0^(-j) * dot(∇f(x), z_k - x) # Armijo linesearch 
       #println(j, " ", stptest)
 
       if stptest > 0.0 # Armijo test   
        j += 1 
        else
        γ = 2.0^(-j)
        break
       end       
    end
    
    β = (- dot(z_k - x, ∇f(x)) * β^2) / (2.0 * (f(x + β * (z_k - x)) - f(x) - β * dot(z_k - x, ∇f(x))))
    
    if β < β1 || β > β2
       β = β / 2.0
    end
    return (γ, β)
 
    println("Step length too small!.")
 end

 function method1(f, ∇f, x, ε, max_iter, GPA1)

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
        newsteps = GPA1(x, f, ∇f, projection, σ, imax_iter, γ_inicial, β_inicial) # GPA1
        z = projection(x - newsteps[2]  * ∇f(x)) # GPA1
        x = x + newsteps[1] * (z - x) # GPA1
        push!(iteration_time, it)
        push!(stepsizes_β, newsteps[2])
        push!(stepsizes_γ, newsteps[1])
        
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