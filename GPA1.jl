## Projected Gradient Method with GPA1 ##

function GPA1(x, f, ∇f, σ, min_step, γ_start, zk)
    ierror = 0
    γ = γ_start
    j = 0
    
    if x == x0
        evalf = 1
        else
        evalf = 0
    end

    while true
     evalf += 1 
     x_plus = x + 2.0^(-j) * (zk - x)   
     stptest = f(x_plus) - f(x) - σ * 2.0^(-j) * dot(∇f(x), zk - x) 
        
        if stptest > 0.0   
           j += 1 
        else
           γ = 2.0^(-j)
           #println("γ = ", γ)
            if γ < min_step
               ierror = 1
               println("Step length too small!")
               break
            end
            break  
        end       
    end
    
    return (γ, ierror, evalf)
end

function method1(x, f, ∇f, ε, max_iter, GPA1, projection)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    evalf_γ = Float64[]
    iteration_time = Float64[]
    
    # Initialization
    ierror = 0
    iter = 0
    seqx = x
    t0 = time()
    β = β2

    if norm(x - projection(x - ∇f(x))) < ε
        # println("x = ", x)
        # println("x - ∇f(x) = ", x - ∇f(x))
        # println("projection(x - ∇f(x)) = ", projection(x - ∇f(x)))
        info = 0
        et = time() - t0
        evalf_γ = 0
        println("x0 is a stationary point!")
        return (x, f(x), info, et, ierror, seqx, evalf_γ)
    end

    while true
        xk = copy(x)
        it0 = time()
        ∇fx = ∇f(x)
        gradnorm = norm(∇fx)

        zk = projection(x - β * ∇f(x))

        β = (-dot(zk - x, ∇f(x)) * β^2) / (2.0 * (f(x + β * (zk - x)) - f(x) - β * dot(zk - x, ∇f(x))))
        
        if β < β1
        β = β1
        elseif β > β2
        β = β2
        end

        (γ, ierror, evalf) = GPA1(x, f, ∇f, σ, min_step, γ_start, zk)
        
        if ierror == 1
            break
        end 

        x = x + γ * (zk - x)
        #println("x = ", x)
        seqx = [seqx x] 
        it = time() - it0
        
        push!(iteration_time, it)
        push!(stepsizes_β, β)
        push!(stepsizes_γ, γ)
        push!(evalf_γ, evalf)
        push!(fvals, f(x)) 
        push!(gradnorms, gradnorm)
        
        # First stopping condition
        if norm(x - xk) < ε
           println("The solution has found!")
           break
        end

        # Update iteration
        iter += 1

        # Second stopping condition
        if iter > max_iter
            println("Maximum of iterations was achieved! Stopping...")
            ierror = 2
            break
        end
    end
    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ, evalf_γ = evalf_γ, iteration_time = iteration_time)
    et = time() - t0
    return (x, f(x), info, et, ierror, seqx, sum(evalf_γ))
end