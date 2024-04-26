# Busca de Armijo ao longo das direções viáveis

using LinearAlgebra

function GPA1(x, f, ∇f, projection, σ, imax_iter, γ_inicial, β_inicial)
   
    β = β_inicial
    γ = γ_inicial
    j = 0
    z_k = projection(x - β * ∇f(x)) # Calcula o ponto z_k
   
   #  println("x = ", " ", x)
   #  println("∇f(x) = ", " ", ∇f(x))
   #  println(x - β * ∇f(x))
   #  println("z_k = ", " ", z_k)
   #  println("z_k - x = ", " ", z_k - x)
   #  println("É direção de descida? --->", " ", dot(∇f(x), z_k - x))

    while j < imax_iter 
     stptest = f(x + 2.0^(-j) * (z_k - x)) - f(x) - σ * 2.0^(-j) * dot(∇f(x), z_k - x) # Testa Armijo para o z_k obtido 
       #println(j, " ", stptest)
 
       if stptest > 0.0 # Se a condição de Armijo não for satisfeita, testa o próximo j   
        j += 1 
        else
        γ = 2.0^(-j)
        break
       end       
    end
    β = (- dot(z_k - x, ∇f(x)) * β^2) / (2.0 * (f(x + β * (z_k - x)) - f(x) - β * dot(z_k - x, ∇f(x))))
    println("β = ",β)
    if β < β1 || β > β2
       β = β / 2.0
       #println("β = ", " ", β)
    end
    return (γ, β)
 
    println("Número máximo de iteradas atingido.")
 end

# x = [1, 1.9] # Ponto inicial
# n = length(x0) # Dimensão da Função Dixon-Price
# y = [1, 1] # Centro do conjunto viável (Bola)
# δ = 1.0 # Raio da Bola
# σ = 1.e-4 # Parâmetro da Busca de Armijo
# β_inicial = 1.0 # Comprimento de passo inicial
# γ_inicial = 1.0 
# β1 = 1.e-6 
# β2 = 1.0
# imax_iter = 100 # Máximo de iteradas da estratégia

# par_ordenado = GPA1(x, f, ∇f, projection, σ, imax_iter, γ_inicial, β_inicial)