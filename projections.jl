#
# Funções de projeções para diferentes conjuntos
#

# ## O CONJUNTO C É UMA BOLA DE RAIO δ E CENTRO y ##

function projection(x)
    norm_yx = norm(y-x)
    if norm_yx <= δ
        return x
    else
        projection = y - (δ*(y-x)/norm_yx)
        return projection
    end
end

# println("Projeção de x em C: ", projection(x))
# δ = 1.0
# y = [1,1]
# k = [0,2]
# println(projection(k))