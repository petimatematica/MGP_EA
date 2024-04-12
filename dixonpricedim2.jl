using LinearAlgebra

function f(x)
    return (x[1] - 1)^2 + 2 * (2 * x[2]^2 - x[1])^2
end

function ∇f(x)
    df_dx = 6*x[1] - 8*x[2]^2 - 2
    df_dy = 32*x[2]^3 - 16*x[1]*x[2]
    return [df_dx, df_dy]
end

function projf(x)
    return [x[1],x[2]]-(dot([x[1],x[2]], [-1/sqrt(2),1]))/(3/2)*[-1/sqrt(2),1]
end

x = [1; 1/sqrt(2)]

# println("f(x) = ", f(x))
# println("∇f(x) = ", ∇f(x))
# println("projf(x) = ", projf(x))
# println("produtoint = ", dot([-1/sqrt(2),1], projf(x)))