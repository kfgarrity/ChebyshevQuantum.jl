

using Optim

#function find_energy(V, L, initial)
function find_energy(z0, ini)

#    a = L/2
#    z0 = a * sqrt(2*V)
    
    function f(z)

        L = sqrt.( abs.((z0./z).^2 .- 1))
        R = tan.(z)
        return sum( (L - R).^2)
        #        alpha = sqrt.( 2.0 / L^2 * abs.(V .+ E))
#        k = sqrt.(abs.(2.0* E))

#        return sum(( k .* tan.(k*L/2.0) .- alpha ).^2)
    end

    ret = optimize(f, ini, BFGS())
    println(ret)
    println()
    println("min val ", f(Optim.minimizer(ret)))
    println(Optim.minimizer(ret)) 

    return Optim.minimizer(ret)
    
end
