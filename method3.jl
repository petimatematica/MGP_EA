## Projected Gradient Method with GPA3 ##

function method3(x0, f, ∇f, ε, max_iter, min_step)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    iteration_time = Float64[]
    
    # Initialization
    ierror = 0
    iter = 1
    x = x0
    seqx = x
    β = β_inicial
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
        k = iter
        α = 1/k   
        β = α / gradnorm

        if β < min_step
            ierror = 1
            println("Step length too small!")
            break
        end
        
        x = projection(x - β * ∇fx)
        seqx = [seqx x] 
        it = time() - it0

        push!(fvals, fx)
        push!(gradnorms, gradnorm)
        push!(stepsizes_β, β)
        push!(iteration_time, it)
        push!(stepsizes_γ, 1.0)
        
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
            ierror = 2
            break 
        end
    end
    avalf = fill(0, size(seqx, 2))
    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ, iteration_time = iteration_time)
    et = time() - t0
    return (x, f(x), info, et, ierror, seqx, avalf)
end