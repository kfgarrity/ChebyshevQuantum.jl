module Chebyshev

using ..Interp:cheb
using ..EdgeDetect:find
using ..Interp:make_cheb
using ..Interp:integrateCC
using ..Interp:derivative

struct Cheb{T}

    a::Float64
    b::Float64
    num_sec::Int64
    ranges::Array{Float64,1}
    clist::Array{cheb{T}}
    
end

Base.show(io::IO,C::Cheb) = begin
    println(io, "Cheb object, num_sec=$(C.num_sec-1), ranges=$(round.(C.ranges,digits=7)), limits $(C.a) to $(C.b)")
end


function (C::Cheb)(x)
    if x < C.a || x > C.b
        println("ERROR, $x outside range $(C.a) to  $(C.b)")
        return 0.0
    end

    for i = 2:C.num_sec
        if x <= C.ranges[i]
            return C.clist[i-1](x)
        end
    end
end

Base.size(c::Cheb) = (c.num_sec-1,)
Base.IndexStyle(::Type{<:Cheb}) = IndexLinear()
Base.similar(a::Cheb, ::Type{T}, d::Dims{1}) where {T} = cheb(T, d[1])
Base.getindex(a::Cheb, i::Int) = a.clist[i]
Base.getindex(a::Cheb, i::Tuple) = a.clist[i]
Base.getindex(a::Cheb, i::UnitRange) = a.clist[i[1]]
Base.setindex!(a::Cheb, v, i) = (a.clist[i] = v)
Base.lastindex(a::Cheb) = a.num_sec-1
Base.length(a::Cheb) = a.num_sec-1

Base.:(==)(c1::Cheb, c2::Cheb) = begin

    if c1.num_sec == c2.num_sec && sum(abs.(c1.ranges - c2.ranges)) < 1e-5

        for (c1a, c2a) in zip(c1.clist, c2.clist)
            if c1a != c2a
                return false
            end
        end
        return true
    end
    return false
    
end


function apply(F::Function, c1::Cheb, n::Number)

    if c1[1].thr < 0
        ct = typeof(c1[1])[]
        for cc in c1.clist
            push!(ct, F(cc, n))
        end
        return make_Cheb(c1.ranges, ct, a=c1.a, b = c1.b)
    else
        c1f = x->c1(x)
        function f(x)
            return F(c1f(x),n)
        end
        return make_Cheb(f, thr=c1[1].thr, a=c1.a, b = c1.b)
    end

end

function apply(F::Function, c1::Cheb)

    if c1[1].thr < 0
        ct = typeof(c1[1])[]
	for cc in c1.clist
            push!(ct, F(cc))
	end
        return make_Cheb(c1.ranges, ct, a=c1.a, b = c1.b)
    else
        
        T = typeof(c1[1])
        clist2 = T[]
        for var in c1.clist
            push!(clist2, F(var))
        end
        make_Cheb(c1.ranges, clist2, a=c1.a, b = c1.b)
    end
end

function apply(F::Function, c1::Cheb, c2::Cheb)
    if !(c1.a ≈ c2.a) || !(c1.b ≈ c2.b)
        println("ERROR max min ranges don't match")
        return c1
    end
    
    if c1[1].thr < 0 && c2[1].thr < 0 && c1.num_sec == c2.num_sec

        ct = typeof(c1[1])[]
        for (cc1,cc2) in zip(c1.clist, c2.clist)
            push!(ct, F(cc1, cc2))
        end
        return make_Cheb(c1.ranges, ct, a=c1.a, b = c1.b)

    else
        
        T = typeof(c1[1])
        clist2 = T[]
        c1f = x->c1(x)
        c2f = x->c2(x)
        
        function f(x)
            return F(c1f(x), c2f(x))
        end
        if c1[1].thr < 0 && c2[1].thr < 0
            return make_Cheb(f, N=max(c1[1].N, c2[1].N), thr=-1, a=c1.a, b = c1.b)
        else
            return make_Cheb(f, thr=max(min(c1[1].thr, c2[1].thr),1e-12), a=c1.a, b = c1.b)
        end

    end
    
end


Base.real(c1::Cheb) = begin
    apply(real, c1)
end

Base.imag(c1::Cheb) = begin
    apply(imag, c1)
end

Base.:^(c1::Cheb, n::Number) = begin
    return apply(^, c1, n)
end

Base.:inv(c1::Cheb) = begin
    return apply(^, c1, -1.0)
end



Base.:-(c1::Cheb) = begin
    T = typeof(c1[1])
    clist2 = T[]
    for var in c1.clist
        push!(clist2, -var)
    end
    return make_Cheb(c1.ranges, clist2, a=c1.a, b = c1.b)

end

Base.:*(c1::Cheb, c2::Cheb) = begin
    return apply(*, c1, c2)
end

Base.:+(c1::Cheb, c2::Cheb) = begin
    return apply(+, c1, c2)
end

Base.:-(c1::Cheb, c2::Cheb) = begin
    return apply(-, c1, c2)
end

Base.:/(c1::Cheb, c2::Cheb) = begin
    return apply(/, c1, c2)
end

Base.:*(c1::Cheb, n::Number) = begin
    return apply(*, c1, n)
end

Base.:+(c1::Cheb, n::Number) = begin
    return apply(+, c1, n)
end

Base.:*(n::Number, c1::Cheb ) = c1*n
Base.:+(n::Number, c1::Cheb ) = c1 + n

Base.:-(n::Number, c1::Cheb ) = begin
    return n + (-c1)
end

Base.:-(c1::Cheb, n::Number ) = -c1 + (n)
Base.:/(c1::Cheb, n::Number  ) = c1 * (1.0/n)
Base.:/(n::Number, c1::Cheb) = n*(c1^-1)



function make_Cheb(ranges::Array{Float64,1}, clist; a=-1.0, b=1.0)
    a=Float64(a)
    b=Float64(b)
    num_sec = length(clist)+1
    return Cheb(a,b,num_sec, ranges,clist)
    
end

function make_Cheb(F::Function; N=-1, thr = 1e-12, a=-1.0, b= 1.0, Nsearch=110)

    
    a=Float64(a)
    b=Float64(b)

    T=typeof( F((a+b)/2.0))
    
    edges, contin = find(F, Nsearch, a=a, b = b)
    num_sec = length(contin)+2

    ranges = zeros(num_sec)

    a1 = a
    clist = cheb{T}[]
    ranges[1] = a1
    
    for i = 1:num_sec-2
        ranges[i+1] = (edges[i,1]+edges[i,2])/2.0
        c = make_cheb(F, N=N, thr=thr, a=a1,b=edges[i,1])
        push!(clist, c)
        a1 = edges[i,2]
    end

    c = make_cheb(F,N=N, thr=thr, a=a1,b=b)
    push!(clist, c)
    ranges[end]=b
    
    return Cheb(a,b,num_sec, ranges, clist)
    
end

function integrateCC(C::Cheb)

    ret = integrateCC(C[1])
    if length(C) > 1
        for c in C[2:end]
            ret += integrateCC(c)
        end
    end
    return ret
end

Base.sum(C::Cheb) = begin
    return integrateCC(C)
end

function derivative(C::Cheb)
    clist = typeof(C[1])[]
    for c in C.clist
        push!(clist, derivative(c))
    end
    return make_Cheb(C.ranges, clist, a=C.a, b= C.b)
end

Base.adjoint(C::Cheb) = derivative(C)


end #end module
