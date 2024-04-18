# Busca de Armijo ao longo das direções viáveis

function GPA1(x, f, ∇f, projection, σ, imax_iter, γ_inicial)
   
    β = β_inicial
    γ = γ_inicial
    j = 0
    z_k = projection(x - β * ∇f(x)) # Calcula o ponto z_k
    println(z_k - x)
 
    while j < imax_iter 
     stptest = f(x - 2.0^(-j) * (z_k - x)) - f(x) + σ * 2.0^(-j) * dot(∇f(x), z_k - x) # Testa Armijo para o z_k obtido
 
     println(j, " ", stptest)
 
       if stptest > 0.0 # Se a condição de Armijo não for satisfeita, testa o próximo j   
        j += 1 
        else
        γ = 2.0^(-j)
        break
       end 
    end
    βnew = (- dot(z_k - x, ∇f(x)) * β_inicial^2) / (2.0 * (f(x + β_inicial * (z_k - x)) - f(x) - β_inicial * dot(z_k - x, ∇f(x))))
    #println("βnew = ",βnew)
    if βnew < β1 || βnew > β2
       βnew = β / 2.0
    end
    return (γ, βnew)
    β = βnew
 
    println("Número máximo de iteradas atingido.")
 end

# x0 = [1, 1.9] # Ponto inicial
# n = length(x) # Dimensão da Função Dixon-Price
# y = [1, 1] # Centro do conjunto viável (Bola)
# δ = 1.0 # Raio da Bola
# σ = 1.e-4 # Parâmetro da Busca de Armijo
# β_inicial = 1.0 # Comprimento de passo inicial
# γ_inicial = 1.0 
# β1 = 1.e-6 
# β2 = 1.0
# imax_iter = 100 # Máximo de iteradas da estratégia

# par_ordenado = GPA1(x0, f, ∇f, projection, σ, imax_iter, γ_inicial)