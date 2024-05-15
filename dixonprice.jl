# Dixon-Price com dimensão n

function f(x::Vector) # x is a vector of real coordinates
    n = length(x) # Vector dimension
    sum = (x[1] - 1)^2 # First part of the function 
    for i in 2:n  # For each element 'x[i]' of the vector 'x'...
        sum += i * (2*x[i]^2 - x[i-1])^2 # Part of the function involving the summation
    end
    return sum # Sum of the first part with the second part
end

function ∇f(x::Vector) # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    g = zeros(n)    # Creates a vector of n coordinates initialized to 0
    g[1] = 2 * (x[1] - 1) - 4 * (2 * x[2]^2 - x[1])  # Update the first coordinate of the vector with this result
    for i in 2:(n-1)
        g[i] = 2 * i * (2 * x[i]^2 - x[i-1]) * (4 * x[i]) - (i+1) * 2 * (2 * x[i+1]^2 - x[i]) # Update coordinates from 2 to n-1 with these results 
    end
    g[n] = n * 2 * (2 * x[n]^2-x[n-1]) * (4 * x[n]) # Update the last coordinate with this result
    return g
end

#global_minimum = [2^(-(2^i-2)/2^i) for i in 1:n];
# val_minimun = f(global_minimum)

# println("f(x) = ", f(x))
# println("∇f(x) =", ∇f(x))
# println("Global Minimum of f(x) = ", global_minimum)
# println("Valor mínimo de f = ", val_minimun)