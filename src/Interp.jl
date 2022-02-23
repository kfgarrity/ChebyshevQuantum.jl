module Interp

using LinearAlgebra
using Polynomials
using Plots



struct cheb{T}

    N::Int64
    a::Float64
    b::Float64
    thr::Float64
    f::Array{T,1}
    p::Function
    
end

function (c::cheb)(x)
    return c.p(forward(x, c.a, c.b))
end

Base.size(c::cheb) = (c.N+1,)
Base.IndexStyle(::Type{<:cheb}) = IndexLinear()
Base.similar(a::cheb, ::Type{T}, d::Dims{1}) where {T} = cheb(T, d[1])
Base.getindex(a::cheb, i::Int) = a.f[i]
Base.getindex(a::cheb, i::Tuple) = a.f[i[1]]
Base.setindex!(a::cheb, v, i) = (a.f[i] = v)
Base.lastindex(a::cheb) = (a.N+1,)
Base.length(a::cheb) = a.N+1

    
Base.show(io::IO,c::cheb) = begin
    println(io, "cheb object, N=$(c.N), thr=$(c.thr), limits $(c.a) to $(c.b)")
end

Base.:(==)(c1::cheb, c2::cheb) = begin
    if c1.N == c2.N && c1.a ≈ c2.a && c1.b ≈ c2.b && sum(abs.(c1.f - c2.f)) < 1e-10
        return true
    else
        return false
    end
end

Base.:*(c1::cheb, c2::cheb) = begin

    if c1.a != c2.a || c1.b != c2.b
        println("ERROR, cheb ranges not the same")
        return c1
    end

    #fixed N calculation
    if (c1.thr < 0 || c2.thr < 0 ) && c1.N == c2.N
        #return make_cheb( x-> c1(x)*c2(x) , N=c1.N, thr=-1.0, a=c1.a, b = c1.b)
        return make_cheb( c1.f .* c2.f , thr=-1.0, a=c1.a, b = c1.b)
        
    else #choose N
        return make_cheb( x-> c1(x)*c2(x) , thr=max(min(c1.thr, c2.thr), 1e-12), a=c1.a, b = c1.b)
    end
end

Base.:+(c1::cheb, c2::cheb) = begin

    if c1.a != c2.a || c1.b != c2.b
        println("ERROR, cheb ranges not the same")
        return c1
    end

    #fixed N calculation
    if (c1.thr < 0 || c2.thr < 0 ) && c1.N == c2.N
        return make_cheb( c1.f +c2.f , thr=-1.0, a=c1.a, b = c1.b)
    
    else #choose N
        return make_cheb( c1.f+c2.f , thr=max(min(c1.thr, c2.thr), 1e-12), a=c1.a, b = c1.b)
    end
end

Base.:-(c1::cheb, c2::cheb) = begin

    if c1.a != c2.a || c1.b != c2.b
        println("ERROR, cheb ranges not the same")
        return c1
    end

    #fixed N calculation
    if (c1.thr < 0 || c2.thr < 0 ) && c1.N == c2.N
        return make_cheb( c1(x).f - c2.f , N=c1.N, thr=-1.0, a=c1.a, b = c1.b)

    else #choose N
        return make_cheb( c1(x).f - c2.f , thr=max(min(c1.thr, c2.thr), 1e-12), a=c1.a, b = c1.b)
    end
end

Base.:/(c1::cheb, c2::cheb) = begin

    if c1.a != c2.a || c1.b != c2.b
        println("ERROR, cheb ranges not the same")
        return c1
    end

    #fixed N calculation
    if (c1.thr < 0 || c2.thr < 0 ) && c1.N == c2.N
        return make_cheb( x-> c1(x)/c2(x) , N=c1.N, thr=-1.0, a=c1.a, b = c1.b)

    else #choose N
        return make_cheb( x-> c1(x)/c2(x) , thr=max(min(c1.thr, c2.thr), 1e-12), a=c1.a, b = c1.b)
    end
end


Base.:^(c1::cheb, n::Number) = begin

    #fixed N calculation
    if (c1.thr < 0 ) 
        return make_cheb( x-> c1(x)^n , N=c1.N, thr=-1.0, a=c1.a, b = c1.b)

    else #choose N
        return make_cheb( x-> c1(x)^n , thr=c1.thr, a=c1.a, b = c1.b)
    end
end

Base.:+(c1::cheb, n::Number) = begin

    #fixed N calculation
    if (c1.thr < 0 ) 
        return make_cheb( x-> c1(x)+n , N=c1.N, thr=-1.0, a=c1.a, b = c1.b)

    else #choose N
        return make_cheb( x-> c1(x)+n , thr=c1.thr, a=c1.a, b = c1.b)
    end
end

Base.:+( n::Number, c1::cheb) = c1 + n


Base.:-(c1::cheb) = begin

    #fixed N calculation
    if (c1.thr < 0 ) 
        return make_cheb( x-> -c1(x) , N=c1.N, thr=-1.0, a=c1.a, b = c1.b)

    else #choose N
        return make_cheb( x-> -c1(x) , thr=c1.thr, a=c1.a, b = c1.b)
    end
end

Base.:-(c1::cheb, n::Number) = begin
    return c1 + (-n)
end

Base.:-( n::Number, c1::cheb) = (-c1)+n

Base.:*(c1::cheb, n::Number) = begin

    return make_cheb(  c1.f*n , thr=c1.thr, a=c1.a, b = c1.b)
    
end

Base.:*( n::Number, c1::cheb) = c1 * n

Base.:/(c1::cheb, n::Number) = begin

    return make_cheb(  c1.f/n , thr=c1.thr, a=c1.a, b = c1.b)

end


Base.:/( n::Number, c1::cheb) = begin
    #fixed N calculation
    if (c1.thr < 0 ) 
        return make_cheb( x-> n * c1(x)^-1.0 , N=c1.N, thr=-1.0, a=c1.a, b = c1.b)
        
    else #choose N
        return make_cheb( x-> n * c1(x)^-1.0 , thr=c1.thr, a=c1.a, b = c1.b)
    end
end




function forward(x,a,b)
    return -1.0 .+ (a .- x)/(a-b)*2.0
end
function rev(x,a,b)
    return a .+ (x .- -1.0)*(b-a)/2.0
end

function make_cheb(F::Function; N = -1, thr = 1e-12, a = -1.0, b = 1.0)

    a=Float64(a)
    b=Float64(b)
    thr=Float64(thr)
    
    if b < a
        temp = b
        b = a
        a = b
    end

    if N > 0
        pts = getCpts(N)
        fpts = F.(rev(pts,a,b))
        p = get_interp(fpts, pts=pts)
        return cheb(N, a, b, thr, fpts, p)
    else
        function F2(x)
            return F(rev(x,a,b))
        end
        p, N = get_interp(F2, thr=thr)
        pts = getCpts(N)
        fpts = F.(rev(pts,a,b))
        return cheb(N, a, b, thr, fpts, p)
        
    end
    
end
    
function make_cheb(fpts::Array; a=-1.0, b= 1.0, thr = 1e-12)

    N = length(fpts) - 1
    p = get_interp(fpts)
    cheb(N, a, b, thr, fpts, p)
    
end
    

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

function derivative(c::cheb, n=1, pts=missing, D=missing)
    N = c.N
    if ismissing(pts)
        pts = getCpts(N)
    end
    if ismissing(D)
        D =  getD(pts)
    end
    fnew = D^n * c.f / ((c.b-c.a)/2.0)^n
    
    return make_cheb(fnew, a=c.a, b=c.b, thr=c.thr)

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

    p = get_interp(f, 1)
    
    for N = vcat(1:50,51:2:205, 210:10:2000)

        p = get_interp(f, N)
        if maximum(abs.(p.(x) - fx)) < thr
            println("approx using $N evals to err $thr")
            return p, N
        end
    end
    println("best found is 2000, err is ", maximum(abs.(p.(x) - fx)))
    return p, 2000
    
end




function getW(N)

    Nd2 = Int64(floor(N/2))
    pts = getCpts(N)
    w = zeros(N+1)
    z = zeros(Nd2*2+1)

    for m = 0:Nd2
        z[:] .= 0.0
        z[2*m+1] = 1.0
        T2m = ChebyshevT(z)

        if m == 0 || m == Nd2
            f2 = 0.5  / (1 - 4*m^2)
        else
            f2 = 1.0  / (1 - 4*m^2)
        end
        
        for j = 1:(N+1)
            
            if j == 1 || j == N+1
                f = 2.0/N
            else
                f = 4.0/N
            end
    
        
            w[j] += f*f2*T2m(pts[j])
            
        end
    end
    return w, pts
end

function integrateCC(f::Function, N::Integer; w = missing, pts = missing, a = -1.0, b = 1.0 )

    if ismissing(w) || ismissing(pts)
        w, pts = getW(N)
    end

    return sum(f.(rev(pts,a,b)) .* w) * (b - a)/2.0

end

function integrateCC(f::Array; w = missing, pts = missing, a = -1.0, b = 1.0 )

    N = length(f)-1
    if ismissing(w) || ismissing(pts)
        w, pts = getW(N)
    end

    return sum(f .* w)*(b - a)/2.0

end

function integrateCC(c::cheb; w=missing, pts=missing)

    #    return integrateCC(c, w=w, pts=pts)*(c.b - c.a)/2.0
    return integrateCC(c.f, w=w, pts=pts, a= c.a, b = c.b)

end

Base.sum(a::cheb) = integrateCC(a)
Base.adjoint(a::cheb) = derivative(a)


#=
function integrate(f::Function, N::Integer)
    
    pts = getCpts(N)
    fp = f.(pts)

#    Nmax = Int64(floor(N/2))
#    println("Nmax $Nmax")
#    println("a ", size(fp[3:2:end]))
#    println("aa ", size(2:2:N))
#    println("b ", size(1.0 ./ (1 .- (2 * (1:Nmax)).^2) ))
#    
#    int = fp[1] 
#    for k = 1: Int(floor(N/2))
#        int += 2 * fp[2*k] / (1 - (2*k)^2)
#    end
    #    return int

    println("pts ", pts)
    
    t = pi *(1:N) / (N+1)

    println("cos t", cos.(t))

    
    f1 = sin.(t) * 2 / (N+1)
    
    w = zeros(N)
    for j = 1:N
        for m = 1:N
            w[j] += sin.(t[j]) * 2 / (N+1) * sin(m*t[j]) * (1 .- cos.(m*pi)) / m
        end
    end
    #    w = f1 * sum(f2)

    int = sum( f.(cos.(t)) .* w)

    return int
    
end
=#

function get_interp(f::Function, N::Integer; pts=missing)
    if ismissing(pts)
        pts = getCpts(N)
    end
    
    fpts = f.(pts)
    return get_interp(fpts,  pts=pts)
    
end

function get_interp(fpts::Array; pts=missing)

    N = length(fpts)-1
    
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
    

function myplot(c::cheb; npts = 200)

    x = (c.a:(c.b-c.a)/npts:c.b)
    
    plot(x, c.(x) )
    
end

function myplot(f::Function, c::cheb; npts = 200)

    myplot(c,f,npts=npts)
    
end

function myplot(c::cheb, f::Function; npts = 200)

    x = (c.a:(c.b-c.a)/npts:c.b)
    
    plot(x, f.(x), label="function", linewidth=4 )
    plot!(x, c.(x), label="cheb", linewidth=2 )
    
end

function myplot_err(f::Function,c::cheb; npts = 200)

    myplot_err(c,f,npts=npts)
end

function myplot_err(c::cheb, f::Function; npts = 200)

    x = (c.a:(c.b-c.a)/npts:c.b)
    
    plot(x, c.(x) -  f.(x) , label="error")

    
end

end #end module
