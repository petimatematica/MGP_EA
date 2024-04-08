using LinearAlgebra

x = [1,2]
n=length(x)

function f(x) 

    return (x[1]-1)^2 + sum(i*((2*x[i]^2-x[i-1]))^2 for i in 2:n)

end

function grad_f()
    
end

global_minimum = [2^(-(2^i-2)/2^i) for i in 1:n]

println("f(x) = ", dixon_price(x))
println("Global Minimum of f(x) = ", global_minimum)