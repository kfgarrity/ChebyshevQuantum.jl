module Op
using LinearAlgebra
using ..Interp:getD
using ..Interp:getCpts
using ..Interp:getCpts1
using ..Interp:baryW1
using ..Interp:barymat
using ..Interp:forward
using ..Interp:rev
using ..EdgeDetect:find

function setup(; N = 20, bc1 = [0,0], bc2 = [0,0], ranges=zeros(0,2), contin = Bool[], A0=0.0, A1=0.0, A2=0.0, B=0.0, a = -1.0, b = 1.0)

    if typeof(A0) <: Number
        function A0f(x); return A0; end
    else
        A0f = A0
    end

    if typeof(A1) <: Number
        function A1f(x); return A1; end
    else
        A1f = A1
    end

    if typeof(A2) <: Number
        function A2f(x); return A2; end
    else
        A2f = A2
    end

    if typeof(B) <: Number
        function Bf(x); return B; end
    else
        Bf = B
    end

    println("asdf")
    println("typeof(A2f) ", typeof(A2f))
    return setup_fn(; N=N, bc1=bc1, bc2=bc2, ranges=ranges, contin = contin, A0=A0f, A1=A1f, A2=A2f, B=Bf, a = a, b = b)
    return 0.0

end

function getmats(N; a=-1.0, b = 1.0)

    D = getD(N);
    pts = getCpts(N);

    Dab = D/ ((b-a)/2.0)
    pts_f = rev(pts,a,b)

    return pts_f, Dab

end


function setup_fn(; N = 20, bc1 = [0,0], bc2 = [0,0], ranges=zeros(0,2), contin = Bool[], A0::Function=0.0, A1::Function=0.0, A2::Function=0.0, B::Function=0.0, a = -1.0, b = 1.0)

    println("b")
    
#    ranges = ranges[contin,:]

    if typeof(N) <: Number
        N = ones(Int64, 1+size(ranges)[1]) * N
    end
    
    DL = []
    VL = []

    a1 = a
    b1 = a
    for c in 1:size(ranges)[1]
        b1 = ranges[c,1]
        pts, D = getmats(N[c]; a = a1+1e-15, b = b1-1e-15)
        push!(DL, D)
        push!(VL, pts)
        a1 = ranges[c,2]
        
    end
    b1 = b
    Vpts, Dab = getmats( N[end]; a = a1+1e-15, b = b1-1e-15)
    push!(DL, Dab)
    push!(VL, Vpts)
    

    num_parts = length(DL)
    hsize = 0
    for (c,(d,v)) in enumerate(zip(DL,VL))
        hsize += size(d)[1]
    end

    H = zeros(hsize, hsize)
    S = zeros(hsize, hsize)

    rhs = zeros(hsize)
    
    
    Bound = hcat(1.0, zeros(1,hsize-1))  #boundary conditions
    Bound = [Bound; hcat(zeros(1,hsize-1), 1.0)]

    rhs[end-1] = bc1[1]
    rhs[end] = bc2[1]    

    
    n1r = 0
    n2r = 0
    n1c = 0
    n2c = 0
    
    for (c,(d,pts)) in enumerate(zip(DL,VL))

        ptsmat = repeat(pts, 1,size(d)[1])
        

#        println("A2 ", A2.(pts))
        D2 = (A2.(ptsmat)) .*(d*d)
        D1 = (A1.(ptsmat)) .*(d)

        #D2 = A2(1.0)*(d*d)
        #D1 = A1(1.0)*(d)
        
        v = diagm(A0.(pts))
        
        
        N = size(d)[1] - 1
        
        x = getCpts(N[c])
        w = baryW1(N[c])
        y = getCpts1(N[c]-2)
        P = barymat(y, x, w)

        H0 = v + D2 + D1
        
        n2r += size(v)[1] - 2
        n2c += size(v)[1] 
        H[(n1r+1):n2r,(n1c+1):n2c] += P*H0
        S[(n1r+1):n2r,(n1c+1):n2c] += P*I(N+1)

        if c < num_parts #add continuity of psi, psi'
            Bound = [hcat(zeros(1,n2c-1), 1, -1, zeros( 1,hsize-(n2c+1))); Bound] #make interface psi continuous
            st = size(DL[2])[2]
            Bound = [vcat(zeros(n1c,1), DL[c][end,:], -DL[c+1][1,:], zeros(hsize-st-n2c ,1))'; Bound] #make deriviates contin
        end

        n2r = n2r
        n1r = n2r
        n2c = n2c
        n1c = n2c
        
    end

    
    sb = size(Bound)[1]
    H[end-sb+1:end, :] = Bound
    S[end-sb+1:end, :] .= 0.0
    
    return H, S, rhs
end







end #end module
