# Busca de Armijo ao longo da fronteira
include("projections.jl")
include("dixonpricedim2.jl")

function GPA2(x, f, ∇f, projection, σ, imax_iter, β_inicial)
   
   β = β_inicial
   j = 0

   while j < imax_iter 
    z_kj = projection(x - β * 2.0^(-j) * ∇f(x)) # Calcula o ponto z_kj
    stptest = f(z_kj) - f(x) + σ * dot(∇f(x), x - z_kj) # Testa Armijo para o z_kj obtido

    #println(j, " ", stptest)

      if stptest > 0.0 # Se a condição de Armijo não for satisfeita, testa o próximo j   
       j += 1 
       else
       β = β_inicial * 2.0^(-j)
       return β
      end
   
   end

   println("Número máximo de iteradas atingido.")
end

newstep = GPA2(x, f, ∇f, projection, σ, imax_iter, β_inicial)
#println("β = ", newstep)

