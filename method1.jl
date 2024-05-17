# ## Projected Gradient Method with GPA1 ##

function GPA1(x, f, ∇f, projection, σ, min_step, γ_inicial, β_inicial)
   
    β = β_inicial
    γ = γ_inicial
    j = 0
    z_k = projection(x - β * ∇f(x))

    while γ > min_step 
     stptest = f(x + 2.0^(-j) * (z_k - x)) - f(x) - σ * 2.0^(-j) * dot(∇f(x), z_k - x)  
 
       if stptest > 0.0   
        j += 1 
        else
        γ = 2.0^(-j)
        break
       end       
    end
    
    β = (- dot(z_k - x, ∇f(x)) * β^2) / (2.0 * (f(x + β * (z_k - x)) - f(x) - β * dot(z_k - x, ∇f(x)))) # Interpolation
    
    if β < β1 || β > β2
       β = β / 2.0
    end
    return (γ, β)
 
    println("Step length too small!")
end

function method1(x0, f, ∇f, ε, max_iter, GPA1)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    iteration_time = Float64[]
    seqx = Vector{Float64}[]
    
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
        #push!(seqx, x_k)
        it0 = time()
        ∇fx = ∇f(x)
        fx = f(x) 
        gradnorm = norm(∇fx)
        newsteps = GPA1(x, f, ∇f, projection, σ, min_step, γ_inicial, β_inicial)
        z = projection(x - newsteps[2]  * ∇f(x)) 
        x = x + newsteps[1] * (z - x)
        it = time() - it0 
        push!(iteration_time, it)
        push!(stepsizes_β, newsteps[2])
        push!(stepsizes_γ, newsteps[1])
        push!(fvals, fx) 
        push!(gradnorms, gradnorm)
        
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