## Dixon-Price Function ##

function f(x::Vector) 
    n = length(x) 
    sum = (x[1] - 1)^2  
    for i in 2:n  
        sum += i * (2*x[i]^2 - x[i-1])^2 
    end
    return sum 
end

function ∇f(x::Vector) 
    n = length(x)   
    g = zeros(n)    
    g[1] = 2 * (x[1] - 1) - 4 * (2 * x[2]^2 - x[1])  
    for i in 2:(n-1)
        g[i] = 2 * i * (2 * x[i]^2 - x[i-1]) * (4 * x[i]) - (i+1) * 2 * (2 * x[i+1]^2 - x[i])  
    end
    g[n] = n * 2 * (2 * x[n]^2-x[n-1]) * (4 * x[n])
    return g
end

#global_minimum = [2^(-(2^i-2)/2^i) for i in 1:n];
# val_minimun = f(global_minimum)
