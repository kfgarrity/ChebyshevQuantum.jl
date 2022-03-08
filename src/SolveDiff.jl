module DiffEq

using LinearAlgebra
using ..Op:setup
using ..EdgeDetect:find
using ..Chebyshev:make_cheb


function solve(;N = 20, bc1 = [:a,0,0.0], bc2 = [:b,0,0.0], ranges=zeros(0,2), contin = Bool[], A0=0.0, A1=0.0, A2=0.0, B=0.0, a = -1.0, b = 1.0, dosplit=true)

    
#    if split && !(typeof(A0) <: Number) && length(ranges) == 0
#        ranges, contin = find(A0, a=a,b=b)
#    elseif split && !(typeof(B) <: Number) && length(ranges) == 0
#        ranges, contin = find(B, a=a,b=b)
#    else
#        ranges = zeros(0,2)
#        contin = Bool[]
#    end

    
    H, S, rhs =  setup(N=N, a=a,b=b,bc1=bc1, bc2=bc2, ranges=ranges, contin=contin, A0=A0, A1=A1,A2=A2,B=B,dosplit=dosplit)
    u = H\rhs

    C = make_cheb(u, a = a, b = b)

    return C

end

    

end #end module
