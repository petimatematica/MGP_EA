# Dixon-Price com dimensão n

x = [1,2]
n=length(x)

function f(x) 

    return (x[1]-1)^2 + sum(i*((2x[i]^2-x[i-1]))^2 for i in 2:n)

end

function gradf(x)
    
    gradf_x1 = 2x[1]-2
    gradf_xi = [8i*x[i]*(2x[i]^2-x[i-1])-(2i-2)*(2x[i+1]^2-x[i]) for i in 2:(n-1)]
    gradf_xn = 8n*x[n]*(2x[n]^2-x[n-1])
    return vcat(gradf_x1, gradf_xi, gradf_xn)
    
end

global_minimum = [2^(-(2^i-2)/2^i) for i in 1:n]
val_minimun = f(global_minimum)




println("f(x) = ", f(x))
println("gradf(x) =", gradf(x))
println("Global Minimum of f(x) = ", global_minimum)
println("Valor mínimo de f = ", val_minimun)