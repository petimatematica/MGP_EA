# Busca de Armijo ao longo da fronteira

using LinearAlgebra, DataFrames

include("dixonpricedim2.jl")
include("projections.jl")

σ = 0.5
max_iter = 100

function GPA2(x, f, ∇f, projection)
   
   β_inicial = 10.0
   β = β_inicial
   j = 0

   while j <= max_iter 
   z_kj = projection(x - β * 2.0^(-j) *  ∇f(x)) # Calcula o ponto z_kj
   stptest = f(z_kj) - f(x) + σ * dot(∇f(x), x - z_kj) # Testa Armijo para o z_kj obtido

   if stptest <= 0.0 # Se a condição de Armijo não for satisfeita, testa o próximo j   
   β = β_inicial * 2.0^(-j)
   return β

   else
   j += 1 
   end

   println("Número máximo de iteradas atingido.")
   return β
   
   end
   
end

newstep = GPA2(x, f, ∇f, projection)
println("β =", newstep)