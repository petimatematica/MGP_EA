## Projected Gradient Method with GPA3 ##

 function method3(f, ∇f, ε, max_iter)

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
        k = iter
        α = 1/k   
        β = α / norm(∇f(x))  
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