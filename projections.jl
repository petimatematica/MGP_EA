#
# Funções de projeções para diferentes conjuntos
#

# ## O CONJUNTO C É UMA BOLA DE RAIO δ E CENTRO y ##

function projection(x)
    norm_yx = norm(y-x)
    if norm_yx <= δ
        return x
    else
        projection = y + δ*(y-x) / norm_yx
        return projection
    end
end

# println("Projeção de x em C: ", projection(x))

k = [2,3]
println(projection(k))