using LinearAlgebra

x = [1,2]
n=length(x)

function dixon_price(x) 

    return (x[1]-1)^2 + sum(i*((2*x[i]^2-x[i-1]))^2 for i in 2:n)

end

println("f(x) = ", dixon_price(x))