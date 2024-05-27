## Projected Gradient Method with GPA1 ##

# include("dixonprice.jl")
# include("projections.jl")

function GPA2(x, f, ∇f, projection, σ, min_step, β_start)
    ierror = 0
    β = β_start
    j = 0
 
    while true
     zkj = projection(x - β * 2.0^(-j) * ∇f(x)) 
     stptest = f(zkj) - f(x) + σ * dot(∇f(x), x - zkj) 
 
        if stptest > 0.0 
         j += 1 
        else
            β = β_start * 2.0^(-j)
            if β < min_step
               ierror = 1
               println("Step length too small!")
               break
            end
            break
        end
    end 
    return (β, ierror, j)
end

# projection = projection1
# x0 = rand(2)
# n = length(x0) 
# σ = 1.e-4 
# ε = 1.e-5 
# β_start = 1.0
# β1 = 1.e-6
# β2 = 1.0
# γ_start = 1.0
# min_step = 1.e-5

# hey = GPA2(x0, f, ∇f, projection, σ, min_step, β_start)

function method2(x0, f, ∇f, ε, max_iter, GPA2)

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
        println("x0 is a stationary point!")
        return (x, f(x))
    end
    
    while true
        xk = copy(x)
        it0 = time()
        ∇fx = ∇f(x) 
        fx = f(x)
        gradnorm = norm(∇fx)
        newstep = GPA2(x, f, ∇f, projection, σ, min_step, β_start)

        if ierror == 1
            break
        end  

        x = projection(x - newstep[1] * ∇fx)
        seqx = [seqx x]
        it = time() - it0
        push!(fvals, fx)
        push!(gradnorms, gradnorm)
        push!(iteration_time, it)
        push!(stepsizes_γ, 1.0)
        push!(stepsizes_β, newstep[1]) 
        push!(avalf, newstep[3])
        
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

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, avalf_β = avalf, stepsizes_γ = stepsizes_γ, iteration_time = iteration_time)
    et = time() - t0
    return (x, f(x), info, et, ierror, seqx, avalf)
end