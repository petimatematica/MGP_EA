## Projected Gradient Method with GPA1 ##

function GPA1(x, f, ∇f, projection, σ, min_step, γ_start, β_start)
    ierror = 0
    β = β_start
    γ = γ_start
    j = 0
    zk = projection(x - β * ∇f(x))

    while true 
     stptest = f(x + 2.0^(-j) * (zk - x)) - f(x) - σ * 2.0^(-j) * dot(∇f(x), zk - x)  
        
        if stptest > 0.0   
           j += 1 
        else
           γ = 2.0^(-j)
            if γ < min_step
               ierror = 1
               #println("Step length too small!")
               break
            end
            break  
        end       
    end
    
    β = (-dot(zk - x, ∇f(x)) * β^2) / (2.0 * (f(x + β * (zk - x)) - f(x) - β * dot(zk - x, ∇f(x))))
    
    if β < β1 || β > β2
       β = β / 2.0
    end
    return (γ, β, ierror, j)
end

function method1(x0, f, ∇f, ε, max_iter, GPA1)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    avalf = Float64[]
    iteration_time = Float64[]
    
    # Initialization
    ierror = 0
    iter = 0
    x = x0
    seqx = x
    t0 = time()

    if norm(x - projection(x - ∇f(x))) < ε
        #println("x0 is a stationary point!")
        return (x, f(x))
    end

    while true
        xk = copy(x)
        it0 = time()
        ∇fx = ∇f(x)
        fx = f(x) 
        gradnorm = norm(∇fx)
        newsteps = GPA1(x, f, ∇f, projection, σ, min_step, γ_start, β_start)
        
        if ierror == 1
            break
        end 

        z = projection(x - newsteps[2] * ∇f(x)) 
        x = x + newsteps[1] * (z - x)
        seqx = [seqx x] 
        it = time() - it0
        
        push!(iteration_time, it)
        push!(stepsizes_β, newsteps[2])
        push!(stepsizes_γ, newsteps[1])
        push!(avalf, newsteps[4])
        push!(fvals, fx) 
        push!(gradnorms, gradnorm)
        
        # First stopping condition
        if norm(x - xk) < ε
           #println("The solution has found!")
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
    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ, avalf_γ = avalf, iteration_time = iteration_time)
    et = time() - t0
    return (x, f(x), info, et, ierror, seqx, avalf)
end