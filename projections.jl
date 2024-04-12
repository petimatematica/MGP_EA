#
# Funções de projeções para diferentes conjuntos
#

using LinearAlgebra

## O CONJUNTO C É UMA BOLA DE RAIO δ E CENTRO y ##

function projection(x)
    norm_yx = norm(y-x)
    if norm_yx <= δ
        return x
    else
        return y + δ*(y-x) / norm_yx
    end
end

# Parâmetros de projeção
x = [5,10]
δ = 2.0
y = [0,1]

# println("Projeção de x em C: ", projection(x))

## O CONJUNTO C É UMA BOLA DE RAIO δ E CENTRO y ##