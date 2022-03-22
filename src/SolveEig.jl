module SolveEig

using LinearAlgebra
using ..Interp:getD
using ..Interp:getCpts
using ..Interp:getCpts1
using ..Interp:baryW1
using ..Interp:barymat
using ..Interp:forward
using ..Interp:rev
using ..EdgeDetect:find
using ..Op:setup

#=
function many_ham(V, N; d2 = -0.5, d1 = 0.0, a = -1.0, b = 1.0, Nsearch=103)

    edges, contin = find(V, Nsearch, a=a, b=b)
#    edges = [0.0 0.0]

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
    
    DL = []
    VL = []
    Vpts, Dab = setup_ham2(V, N; a = a, b = -1.0 - 1e-30)    
    push!(DL, Dab)
    push!(VL, Vpts)
    Vpts, Dab = setup_ham2(V, N; a = -1.0+1e-30, b = 1.0-1e-30)    
    push!(DL, Dab)
    push!(VL, Vpts)
    Vpts, Dab = setup_ham2(V, N; a =  1.0+1e-30, b = b)    
    push!(DL, Dab)
    push!(VL, Vpts)

    
#    println("SIZE HL ", size(HL))
    n = 0
    for h in DL
        println( size(h)[1])
        n += size(h)[1]
    end

    
    nh = n - length(DL) + 1 
    #nh = n
    H = zeros(nh, nh)

    Dbig = zeros(nh, nh)
    Vbig = zeros(nh, nh)
    n1 = 0
    n2 = 0

    num_secs = length(DL)
    
    for (c,(d,v)) in enumerate(zip(DL, VL))
        n2 += size(d)[1]
        Dbig[(n1+1):n2,(n1+1):n2] += d
        Vbig[(n1+1):n2,(n1+1):n2] += v
#        println(n1+1, " ", n2)
        n2 = n2-1
        n1 = n2
        
    end

    n1 = 0
    n2 = 0
    for (c,(d,v)) in enumerate(zip(DL, VL))
#        println("c $c")
#        println(v)
        
        n2 += size(d)[1]
        if c != 1
            Vbig[n1+1,n1+1] += v[1,1] / 2.0
        end
        if c != num_secs
            Vbig[n2,n2] = v[end,end] / 2.0
        end
        n2 = n2-1
        n1 = n2
    end

    
    
    D2 = Dbig*Dbig 

    H = d2 * D2 + Vbig
#    H= Vbig
    
#=    n1 = 0
    n2 = 0
    for h in HL
        n2 += size(h)[1]
        H[(n1+1):n2,(n1+1):n2] += h
        println(n1+1, " ", n2)
        n2 = n2
        n1 = n2
        
    end
=#
    
#=    n1 = 0
    n2 = 0
    for h in HL
        n2 += size(h)[1]
        H[(n1+1):n2,(n1+1):n2] += h
        println(n1+1, " ", n2)
        n2 = n2-1
        n1 = n2
        
    end
=#    
    
    return H
end
    
function total_ham_piecewise(V, N, ranges; d2 = -0.5, d1 = 0.0)
    
    DL = []
    VL = []
    for c in 1:length(ranges)-1
        a = ranges[c]
        b = ranges[c+1]
        Vpts, Dab = setup_ham2(V, N; a = a+1e-15, b = b-1e-15)
        push!(DL, Dab)
        push!(VL, Vpts)
    end


    num_parts = length(DL)
    hsize = 0
    for (c,(d,v)) in enumerate(zip(DL,VL))
        hsize += size(d)[1]
    end

    H = zeros(hsize, hsize)
    S = zeros(hsize, hsize)

    B = hcat(1.0, zeros(1,hsize-1))  #boundary conditions
    B = [B; hcat(zeros(1,hsize-1), 1.0)]

    n1r = 0
    n2r = 0
    n1c = 0
    n2c = 0
    
    for (c,(d,v)) in enumerate(zip(DL,VL))
        D2 = d2*(d*d)
        
        N = size(d)[1] - 1
        
        x = getCpts(N)
        w = baryW1(N)
        y = getCpts1(N-2)
        P = barymat(y, x, w)

        H0 = v + D2
        
        n2r += size(v)[1] - 2
        n2c += size(v)[1] 
        
 #       println("a ", size(H[(n1r+1):n2r,(n1c+1):n2c]))
 #       println("p ", size(P))
 #       println("H0 ", size(H0))
        H[(n1r+1):n2r,(n1c+1):n2c] += P*H0
        S[(n1r+1):n2r,(n1c+1):n2c] += P*I(N+1)

        if c < num_parts #add continuity of psi, psi'
#            println("B ", size(B))
#            println("hcat ", size(hcat(zeros(1,n2c-1), 1, -1, zeros(1, hsize-(n2c+1)))))
            B = [B; hcat(zeros(1,n2c-1), 1, -1, zeros( 1,hsize-(n2c+1)))] #make interface psi continuous
            st = size(DL[2])[2]
            B = [B; vcat(zeros(n1c,1), DL[c][end,:], -DL[c+1][1,:], zeros(hsize-st-n2c ,1))'] #make deriviates contin
            #            B = [B; hcat(zeros(1,n1c), d[end,:], -DL[c+1][1,:], zeros(1, hsize - (n2c+1))) ]
        end

        n2r = n2r
        n1r = n2r
        n2c = n2c
        n1c = n2c
        
    end

    
    sb = size(B)[1]
    H[end-sb+1:end, :] = B
    S[end-sb+1:end, :] .= 0.0
    

#    println("size Vpts ", size(Vpts))
#    println("size Dab ", size(Dab))
#    println("size x ", size(x))
#    println("size y ", size(y))
#    println("size P ", size(P))
    

#    println("size PH ", size(PH))
    
#    B = [ hcat(1.0, zeros(1,N));  hcat(zeros(1,N), 1.0)]

#    println("size B ", size(B))
    
#    Hbig = [PH; B]
#    s = size(H0)[1]
#    S = collect(I(size(H0)[1]))
#    PS = P * S
#    Sbig = [PS; zeros(2, s)]

    
    return H, S

    
    
    
    #    Vpts, Dab = setup_ham2(V, N; a = a, b = -1.0 - 1e-30)    
#    push!(DL, Dab)
#    push!(VL, Vpts)
#    Vpts, Dab = setup_ham2(V, N; a = -1.0+1e-30, b = 1.0-1e-30)    
#    push!(DL, Dab)
#    push!(VL, Vpts)

    
    
    

end
=#

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

#=
function setup_ham2(V, N; d2 = -0.5, d1 = 0.0, a = -1.0, b = 1.0)

    D = getD(N);
    pts = getCpts(N);

    Dab = D/ ((b-a)/2.0)
#    Dab2 = Dab*Dab

    pts_f = rev(pts,a,b)

    Vpts = diagm(V.(pts_f))
#    H = d2*Dab2 + d1*Dab + Vpts

    return Vpts, Dab

end

=#

function finite_diff_ham(V, N; d2 = -0.5, d1 = 0.0, a = -1.0, b = 1.0)

    h = (b-a)/N
    grid = a:h:(b+1e-15)

    Vpts = diagm(V.(grid))

    D2 = zeros(length(grid), length(grid))
    for i = 2:(length(grid)-1)
        D2[i,i] = 2.0/h^2*d2
        D2[i,i+1] = -1.0/h^2*d2
        D2[i,i-1] = -1.0/h^2*d2
    end
    D2[1,1] = 1.0/h^2*d2
    D2[1,2] = -1.0/h^2*d2
    D2[end,end] = 1.0/h^2*d2
    D2[end, end-1] = -1.0/h^2*d2

    return Vpts - D2 
    

end

function solve(;N = 30, bc1 =  [:a,0,0.0], bc2 = [:b,0,0.0], ranges=zeros(0,2), contin = Bool[], A0=0.0, A1=0.0, A2=-0.5, B=1.0, a = -1.0, b = 1.0, V = missing, dosplit=true)

    if !ismissing(V)
        A0 = V
    end
       

    H, S, rhs =  setup(N=N, a=a,b=b,bc1=bc1, bc2=bc2, ranges=ranges, contin=contin, A0=A0, A1=A1,A2=A2,B=B,dosplit=dosplit)

    vals, vects = eigen( H, S)

    return vals, vects
    
    #vals, vects = eigen( H[2:end-1,2:end-1]);
    #vals, vects = eigen( H);

#    vects_n = normalize_vects(vects)

    #return vals, vects

end


end #end module
