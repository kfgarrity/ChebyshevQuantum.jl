

function goinv(;ω = 1.0, N = 60, a = -10.0, b = 10.0)

    D = ChebyshevQuantum.Interp.getD(N);
    Dab = D/ ((b-a)/2.0)

    KE = (-0.5*Dab*Dab )[2:N,2:N];

    
    function V(x)
        return 0.5 * ω^2 * x^2
    end

    V_m = ChebyshevQuantum.Interp.make_cheb(V; N = N, a = a, b = b)[2:N]

    vals, vectsR = eigen(KE + diagm(V_m))
    
    #    vals, vectsR = ChebyshevQuantum.SolveEig.solve(V=V, N=N, a= a, b= b, dosplit=false)

    println("target vals " , vals[1:4])
    
    rho_target =  ( real.(vectsR[:,1] .* conj(vectsR[:,1])) )

    β = 1.0
    function V_start(x)
        return 0.5 * β^2 * x^2
    end


    #pts = ChebyshevQuantum.Interp.getCpts(N);

    V_work = ChebyshevQuantum.Interp.make_cheb(V_start; N = N, a = a, b = b)[2:N]

    #V_work = ChebyshevQuantum.Cheb.make_Cheb(V_start; N = 60, a = a, b = b)[2:N]

    H = KE + diagm(V_work)
    
    vals, vects = eigen(H)

    rho = real.(vects[:,1] .* conj(vects[:,1]))

    println("rho diff ", sum(abs.(rho - rho_target)))
    
    println("calc vals " , vals[1:4])

    MAT = zeros(Complex{Float64}, N, N)

    MAT[1:(N-1), 1:(N-1) ] = H - I(N-1)*vals[1]
    MAT[1:(N-1), N] = vects[:,1]
    MAT[N,1:(N-1)] = vects[:,1]'
    
    B = zeros(Complex{Float64}, N)

    println( size(rho_target), " " , size(rho), "  ", size(vects[:,1]))
    B[1:(N-1)] = 4*(rho_target - rho) .* vects[:,1]

    v =  MAT \ B

    return v, rho, rho_target, vects, vectsR
    
    #    return vals

end
