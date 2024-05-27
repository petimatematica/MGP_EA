using LinearAlgebra

#
# PROJECTION FUNCTIONS FOR DIFERENT FEASIBLE SETS 
#

## PROJECTION 01: O CONJUNTO C É UMA BOLA DE RAIO δ E CENTRO y ##

function projection1(x)
    y = fill(0, n) # Centro do conjunto viável (Bola)
    δ = 0.5 # Raio da Bola
    norm_yx = norm(y-x)
    if norm_yx <= δ
        return x
    else
        projection = y - (δ*(y-x)/norm_yx)
        return projection
    end
end

## PROJECTION 02: C É O R^n+ ##

function projection2(x)
    return max.(0, x)
end

## PROJECTION 03: C É UM HIPERCUBO DE DIMENSÃO n ##
## Para cada dimensão, a amplitude é o intervalo [-1, 1]

function projection3(x)
    n = length(x)
    a = fill(-1, n)
    b = fill(1, n)
    projection = similar(x)
    for j in 1:n   
        if x[j] < a[j]
           projection[j] = a[j]
        elseif a[j] <= x[j] <= b[j]
           projection[j] = x[j]
        elseif x[j] > b[j]
            projection[j] = b[j] 
        end
    end
    return projection
end

## PROJECTION 04: C É O CONJUNTO DOS PONTOS QUE SATISFAZEM UMA IGUALDADE COM UM PRODUTO INTERNO DEFINIDO ##

function projection4(x)
    n = length(x)
    a = fill(1, n)
    b = 0
    projection = x + ((b - dot(a, x))/norm(a)^2)*a
    return projection
end

## PROJECTION 05: C É O CONJUNTO DOS PONTOS QUE SATISFAZEM UMA IGUALDADE COM UM PRODUTO INTERNO DEFINIDO ##

function projection5(x)
    n = length(x)
    a = fill(1, n)
    b = 0
    projection = x + ((min(0, b - dot(a, x)))/norm(a)^2)*a
    return projection
end

## PROJECTION 06: C É O CONJUNTO DOS PONTOS QUE SATISFAZEM UMA IGUALDADE COM UM PRODUTO DE VETOR E MATRIZ ##

function matrix(rank, n)
    matrix_0 = rand(rank, n)  
    U, Σ, V = svd(matrix_0)
    Σ[rank + 1:end] .= 0      
    A = U[:, 1:rank] * Diagonal(Σ[1:rank]) * V[:, 1:rank]'
    return A
end

function projection6(x)
    n = length(x)
    rank = 2  

    global A_dict
    if isempty(A_dict)
        A_dict = Dict{Int, Matrix{Float64}}()
    end

    if !haskey(A_dict, n)
        A_dict[n] = matrix(rank, n)  
    end

    A = A_dict[n]
    b = fill(0, size(A, 1))
    projection = x - A' * inv(A * A') * (A * x - b)
    return projection
end

## TESTE ##

# y = [0.5, 2, 3, 1];
# println("Projeção de x em C: ", projection6(y))
#hey = projection6(y)