module Interp

using LinearAlgebra

function getD(N::Integer)


    if N == 0
        return zeros(1,1)
    end
    
    pts = getCpts(N)
    return getD(pts)
    
end

function getD(pts)

    N = length(pts)-1
    if N <= 0
        return zeros(1,1)
    end

    c = vcat([2.0], ones(N-1), [2.0]) .* (-1).^(0:N)
    
    X = repeat(-pts, 1,N+1)

    dX = X - X'
    
    D = (c * (1.0 ./c)' ./ (dX + I(N+1)))
    
    D = -D + diagm(sum(D, dims=2)[:])

    return D
    
end

function get_derivative(f::Function, N::Integer; n=1, pts=missing, D= missing)
    if ismissing(pts)
        pts = getCpts(N)
    end
    if ismissing(D)
        D =  getD(pts)
    end

    return get_derivative(f.(pts), N, n=n, pts=pts, D=D)
end

function get_derivative(fp::Array, N::Integer; n=1, pts=missing, D= missing)
    
    fpp = (D^n)*fp

    return get_interp(fpp, N, pts=pts)
    
end


function getCpts(N::Integer)

    if N <= 0
        return 0.0
    else
        return -cos.( (0:N)*pi / N)
    end
        
end

function get_interp(f; thr = 1e-12)
    x = -1:.01:1
    fx = f.(x)

    p = get_interp(f, 0)
    
    for N = vcat(1:200, 210:10:2000)

        p = get_interp(f, N)
        if maximum(abs.(p.(x) - fx)) < thr
            println("approx using $N evals to err $thr")
            return p
        end
    end
    println("best found is 2000, err is ", maximum(abs.(p.(x) - fx)))
    return p
    
end
    

function get_interp(f::Function, N::Integer; pts=missing)
    if ismissing(pts)
        pts = getCpts(N)
    end
    
    fpts = f.(pts)
    return get_interp(fpts, N, pts=pts)
    
end

function get_interp(fpts::Array, N::Integer; pts=missing)

    if ismissing(pts)
        pts = getCpts(N)
    end
    
    g = zeros(N+1)

    function gg(x)

        g[:] .= (-1).^(0:N) ./ (x .- pts)

        #        for k = 0:N
        #            g[k+1] = (-1)^k./(x - pts[k+1])
        #        end
        
        g[1] = g[1] * 0.5
        g[N+1] = g[N+1] * 0.5
        return g
    end
    
    function p(x)

        i = argmin(abs.(x .- pts))
        if abs(x - pts[i]) < 1e-10 #if exact match return function value, avoid divde by zero
            return fpts[i]
        else

            g = gg(x)
            #barycentric
            return sum( fpts.*g) / sum(g)

        end
    end

    return p
    
    
end
    


end #end module
