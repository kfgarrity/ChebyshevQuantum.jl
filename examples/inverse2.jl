using LinearAlgebra
using Optim

function goinv2(;ω = 1.0, N = 60, a = -3.0, b = 3.0)

    D = ChebyshevQuantum.Interp.getD(N);
    Dab = D/ ((b-a)/2.0)

    KE = (-0.5*Dab*Dab )[2:N,2:N];

    xvals = ChebyshevQuantum.Interp.rev(ChebyshevQuantum.Interp.getCpts(N), a, b)
    
    
    
    function V(x)
        return 0.5 * ω^2 * x^2  + 0.1 * x^3
    end

    V_m = ChebyshevQuantum.Interp.make_cheb(V; N = N, a = a, b = b)[2:N]

    vals, vects = eigen(KE + diagm(V_m))
    
    #    vals, vectsR = ChebyshevQuantum.SolveEig.solve(V=V, N=N, a= a, b= b, dosplit=false)

    println("target vals " , vals[1:4])
    
    rho_target = sum( ( real.(vects[:,1:2] .* conj(vects[:,1:2])) ), dims=2)

    #weights
    #w = ( 10^6 .+ 10^8 * xvals.^4)[2:N]
    w = 10^7 ./ rho_target

    
    β = 1.0
    function V_start(x)
        return 0.5 * β^2 * x^2
    end


    #pts = ChebyshevQuantum.Interp.getCpts(N);

    V_work = ChebyshevQuantum.Interp.make_cheb(V_start; N = N, a = a, b = b)[2:N]

    #V_work = ChebyshevQuantum.Cheb.make_Cheb(V_start; N = 60, a = a, b = b)[2:N]

    H = zeros(N-1, N-1)
    
    function f(x)
    
        H[:,:] = KE + diagm(x)
        vals, vects = eigen(H)
        rho = sum(real.(vects[:,1:2] .* conj(vects[:,1:2])), dims=2)
        println("rho diff ", sum(abs.(rho - rho_target)))
#        println("calc vals " , vals[1:4])

        return   sum(w.*(rho - rho_target).^2)
    end

    grad = zeros(N-1)

    MAT = zeros(Complex{Float64}, N, N)
    B = zeros(Complex{Float64}, N)
    rho = zeros(Float64, N)
    
    function g(grad, x)

        H[:,:] = KE + diagm(x)
        vals, vects = eigen(H)

        rho[:] = sum(real.(vects[:,1:2] .* conj(vects[:,1:2])), dims=2)

        grad .= 0.0
        for ii = 1:2

            MAT[1:(N-1), 1:(N-1) ] = H' - I(N-1)*vals[ii]
            MAT[1:(N-1), N] = 2.0*vects[:,ii]
            MAT[N,1:(N-1)] = vects[:,ii]'
        

            #        println( size(rho_target), " " , size(rho), "  ", size(vects[:,1]))
            B[1:(N-1)] = w .* 4.0 .*(rho_target - rho) .* vects[:,ii]

            v =  MAT \ B

            grad[:] += real.(v[1:N-1] .* vects[:,ii])
        end
#        println("grad ", sum(abs.(grad)))
        
        return grad
    end

    ff = f(V_work)
    gg = g(grad, V_work)
    
    println("f ", ff)
    println("g ", gg)

#    return  ff, V_work, V_m, f, g

    opts = Optim.Options( f_tol = 1e-10, g_tol = 1e-10, iterations = 1000, store_trace = true, show_trace = false)    
    ret = optimize(f, g, V_work, BFGS(), opts)
#    ret = optimize(f, V_work, BFGS())

    println("ret")
    println(ret)
    println("min")
    themin = Optim.minimizer(ret)
    
    return ret, themin, V_m, f, g, rho_target, rho
#    return v, rho, rho_target, vects, vectsR
    
    #    return vals

end
