## Projected Gradient Method with GPA1 ##

function GPA2(x, f, ∇f, projection, σ, min_step, β_start)
    ierror = 0
    β = β_start
    j = 0
    
    if x == x0
        evalf = 1
        else
        evalf = 0
    end
    
    zkj = 0

    while true
     evalf += 1 
     zkj = projection(x - β * 2.0^(-j) * ∇f(x))
     stptest = f(zkj) - f(x) + σ * dot(∇f(x), x - zkj) 
 
        if stptest > 0.0 
         j += 1 
        else
            β = β_start * 2.0^(-j)
            #println(β)
            if β < min_step
               ierror = 1
               println("Step length too small!")
               break
            end
            break
        end
    end 
    return (β, zkj, ierror, evalf)
end

function method2(x, f, ∇f, ε, max_iter, GPA2, projection)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    evalf_β = Float64[]
    iteration_time = Float64[]
    
    # Initialization
    ierror = 0
    iter = 0
    seqx = x
    t0 = time()

    if norm(x - projection(x - ∇f(x))) < ε
        info = 0
        et = time() - t0
        evalf_β = 0
        println("x0 is a stationary point!")
        return (x, f(x), info, et, ierror, seqx, evalf_β)
    end
    
    while true
        xk = copy(x)
        it0 = time()
        ∇fx = ∇f(x) 
        fx = f(x)
        gradnorm = norm(∇fx)
        (β, zkj, ierror, evalf) = GPA2(x, f, ∇f, projection, σ, min_step, β_start)

        if ierror == 1
            break
        end  

        x = zkj
        #println("x =", x)
        seqx = [seqx x]
        it = time() - it0
        push!(fvals, fx)
        push!(gradnorms, gradnorm)
        push!(iteration_time, it)
        push!(stepsizes_β, β) 
        push!(evalf_β, evalf)

        #println("Critério de parada -> ", norm(x - xk))
        
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

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, evalf_β = evalf_β, iteration_time = iteration_time)
    et = time() - t0
    return (x, f(x), info, et, ierror, seqx, sum(evalf_β))
end