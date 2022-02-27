module SolveEig

using LinearAlgebra
using ..Interp:getD
using ..Interp:getCpts
using ..Interp:forward
using ..Interp:rev
using ..EdgeDetect:find

function many_ham(V, N; d2 = -0.5, d1 = 0.0, a = -1.0, b = 1.0, Nsearch=103)

#    edges, contin = find(V, Nsearch, a=a, b=b)
    edges = [0.0 0.0]

    a1=a
    HL = []
               
    for i in 1:size(edges)[1]
        b1 = edges[i,1]
        h = setup_ham(V, N, a=a1, b = b1, d2=d2, d1=d1)
        push!(HL, h)
        a1=edges[i,2]
    end
    h = setup_ham(V, N, a=a1, b = b, d2=d2, d1=d1)
    push!(HL, h)
    
    return HL
#    for i in range(edges
#    H1 = setup_ham(V, a=a
    

end


function total_ham(V, N; d2 = -0.5, d1 = 0.0, a = -1.0, b = 1.0, Nsearch=103)

    HL = []
    push!(HL, setup_ham(V, N; a = -1.0, b = 0.0))
    push!(HL, setup_ham(V, N; a = 0.0, b = 1.0))
    #    HL = many_ham(V, N; d2 = d2, d1 = d1, a = a, b = b, Nsearch=Nsearch)

    n = 0
    for h in HL
        println( size(h)[1])
        n += size(h)[1]
    end

    
    nh = n - length(HL) + 1 
    #nh = n
    H = zeros(nh, nh)

#    n1 = 0
#    n2 = 0
#    for h in HL
#        n2 += size(h)[1]
#        H[(n1+1):n2,(n1+1):n2] += h
#        println(n1+1, " ", n2)
#        n2 = n2
#        n1 = n2
#        
#    end

    
    n1 = 0
    n2 = 0
    for h in HL
        n2 += size(h)[1]
        H[(n1+1):n2,(n1+1):n2] += h
        println(n1+1, " ", n2)
        n2 = n2-1
        n1 = n2
        
    end
    
    return H
end
    

    
function setup_ham(V, N; d2 = -0.5, d1 = 0.0, a = -1.0, b = 1.0)

    D = getD(N);
    pts = getCpts(N);

    Dab = D/ ((b-a)/2.0)
    Dab2 = Dab*Dab

    pts_f = rev(pts,a,b)

    Vpts = diagm(V.(pts_f))
    H = d2*Dab2 + d1*Dab + Vpts

    return H

end

function solve(H)

    vals, vects = eigen( H[2:end-1,2:end-1]);
    #vals, vects = eigen( H);

#    vects_n = normalize_vects(vects)

    return vals, vects

end


end #end module
