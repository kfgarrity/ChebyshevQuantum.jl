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

function setup(; N = 20, bc1 = [:a,0,0.0], bc2 = [:b,0,0.0], ranges=zeros(0,2), contin = Bool[], A0=0.0, A1=0.0, A2=0.0, B=0.0, a = -1.0, b = 1.0, dosplit=true)

    a = Float64(a)
    b = Float64(b)
    
    if typeof(A0) <: Number
        function A0f(x); return A0; end
        ranges0 = zeros(0,2)
        contin0 = Bool[]
    else
        A0f = A0
        if dosplit && length(ranges) == 0
            ranges0, contin0 = find(A0f,a=a,b=b)
        end
    end

    if typeof(A1) <: Number
        function A1f(x); return A1; end
        ranges1 = zeros(0,2)
        contin1 = Bool[]
    else
        A1f = A1
        if dosplit && length(ranges) == 0
            ranges1, contin1 = find(A1f,a=a,b=b)
        end
    end

    secondorder=true
    if typeof(A2) <: Number
        if abs(A2) < 1e-15
            secondorder=false
        end
        function A2f(x); return A2; end
        ranges2 = zeros(0,2)
        contin2 = Bool[]
    else
        A2f = A2
        if dosplit && length(ranges) == 0
            ranges2, contin2 = find(A2f,a=a,b=b)
        end
    end

    if typeof(B) <: Number
        function Bf(x); return B; end
        rangesB = zeros(0,2)
        continB = Bool[]
    else
        Bf = B
        if dosplit && length(ranges) == 0
            rangesB, continB = find(Bf,a=a,b=b)
        end
    end

    if length(ranges) == 0 && dosplit

        if dosplit && length(ranges) == 0
            #merge ranges
            ranges = [ranges0; ranges1;ranges2;rangesB]
            contin = [contin0; contin1; contin2; continB]
            s = sortperm(ranges[:,1])
            ranges = ranges[s,:]
            contin = contin[s,:]

#            println("ranges A")
#            println(ranges)
            
            if size(ranges)[1] > 0
                keep = [1]
                for i in 2:(size(ranges)[1])
                    if sum(abs.(ranges[i-1,:] - ranges[i,:])) > 1e-4
                        push!(keep, i)
                    end
                end
                ranges = ranges[keep,:]
                contin = contin[keep,:]
                if length(ranges) > 0
                    println("detected edges")
                    println(ranges)
#                    println("contin")
#                    println(contin)
#                    println()
                end
            end
        end
    end
    
#    println("asdf")
#    println("typeof(A2f) ", typeof(A2f))

    
    
    if secondorder
        return setup_fn(; N=N, bc1=bc1, bc2=bc2, ranges=ranges, contin = contin, A0=A0f, A1=A1f, A2=A2f, B=Bf, a = a, b = b)
    else
        return setup_fn_firstorder(; N=N, bc1=bc1, ranges=ranges, contin = contin, A0=A0f, A1=A1f, B=Bf, a = a, b = b)
    end        

end

function getmats(N; a=-1.0, b = 1.0)

    D = getD(N);
    pts = getCpts(N);

    Dab = D/ ((b-a)/2.0)
    pts_f = rev(pts,a,b)

    return pts_f, Dab

end


function setup_fn(; N = 20, bc1 = [:a,0,0.0], bc2 = [:b,0,0.0], ranges=zeros(0,2), contin = Bool[], A0::Function, A1::Function, A2::Function, B::Function, a = -1.0, b = 1.0)

#    println("b")
    
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

    
    Bound = zeros(0,hsize)

    function apply_bc(bc)
        if ismissing(bc)
            return hcat(1.0, zeros(1,hsize-1)), 0.0
        end            
        if bc[1] == :a && bc[2] == 0
            return hcat(1.0, zeros(1,hsize-1)), bc[3]
        elseif bc[1] == :b && bc[2] == 0
            return hcat(zeros(1,hsize-1), 1.0), bc[3]
        elseif bc[1] == :a && bc[2] == 1
            return hcat(DL[1][1,:]', zeros(1,hsize-size(DL[1])[1])), bc[3]
        elseif bc[1] == :b && bc[2] == 1
#            println("s ", size(zeros(1,hsize-size(DL[end])[1])))
#            println("s2 ",  size(DL[end][end,:]))
            return hcat( zeros(1,hsize-size(DL[end])[1]), DL[end][end,:]'),bc[3]
        else
            println("WARNING, bc not recognized $bc")
            return hcat(1.0, zeros(1,hsize-1)), bc[3]
        end
    end

#    println("bc1 $bc1 bc2 $bc2")
    
    Bound = zeros(0,hsize)
    BC1, v1 = apply_bc(bc1)
    BC2, v2 = apply_bc(bc2)

#    println("hsize $hsize")
#    println("size BC1 ", size(BC1))
#    println("size BC2 ", size(BC2))
    
    Bound = [Bound; BC1]

    if !ismissing(bc1)
        Bound = [Bound; BC2]
    end

    if ismissing(bc2)
        rhs[end] = v1
    else
        rhs[end-1] = v1
        rhs[end] = v2    
    end
    
    
#    Bound = hcat(1.0, zeros(1,hsize-1))  #boundary conditions
#    Bound = [Bound; hcat(zeros(1,hsize-1), 1.0)]




    
    n1r = 0
    n2r = 0
    n1c = 0
    n2c = 0
    
    for (c,(d,pts)) in enumerate(zip(DL,VL))


        
        #        ptsmat = repeat(pts', size(d)[1], 1)

        ptsmat = repeat(pts,1, size(d)[1])
        

#        println("A2 ", A2.(pts))
        D2 = (A2.(ptsmat)) .*(d*d)
        D1 = (A1.(ptsmat)) .*(d)

        #D2 = A2(1.0)*(d*d)
        #D1 = A1(1.0)*(d)
        
        v = diagm(A0.(pts))
        #v = A0.(ptsmat)
        
        
        N = size(d)[1] - 1
        
        x = getCpts(N)
        w = baryW1(N)
        y = getCpts1(N-2)
        P = barymat(y, x, w)

        H0 = v + D2 + D1
        
        n2r += size(v)[1] - 2
        n2c += size(v)[1] 
        H[(n1r+1):n2r,(n1c+1):n2c] += P*H0
        S[(n1r+1):n2r,(n1c+1):n2c] += P*I(N+1)

        rhs[(n1r+1):n2r] = P*B.(pts)
        
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


function setup_fn_firstorder(; N = 20, bc1 = [:a,0,0.0], ranges=zeros(0,2), contin = Bool[], A0::Function, A1::Function, B::Function, a = -1.0, b = 1.0)

    println("solve first order")
    
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

    
    Bound = zeros(0,hsize)

    function apply_bc(bc)
        if ismissing(bc)
            return hcat(1.0, zeros(1,hsize-1)), 0.0
        end            
        if bc[1] == :a && bc[2] == 0
            return hcat(1.0, zeros(1,hsize-1)), bc[3]
        elseif bc[1] == :b && bc[2] == 0
            return hcat(zeros(1,hsize-1), 1.0), bc[3]
        elseif bc[1] == :a && bc[2] == 1
            return hcat(DL[1][1,:]', zeros(1,hsize-size(DL[1])[1])), bc[3]
        elseif bc[1] == :b && bc[2] == 1
#            println("s ", size(zeros(1,hsize-size(DL[end])[1])))
#            println("s2 ",  size(DL[end][end,:]))
            return hcat( zeros(1,hsize-size(DL[end])[1]), DL[end][end,:]'),bc[3]
        else
            println("WARNING, bc not recognized")
            return hcat(1.0, zeros(1,hsize-1)), bc[3]
        end
    end

    println("bc1 $bc1 ")
    
    Bound = zeros(0,hsize)
    BC1, v1 = apply_bc(bc1)

#    println("hsize $hsize")
#    println("size BC1 ", size(BC1))
#    println("size BC2 ", size(BC2))
    
    Bound = [Bound; BC1]

    rhs[end] = v1    

    
    
    n1r = 0
    n2r = 0
    n1c = 0
    n2c = 0
    
    for (c,(d,pts)) in enumerate(zip(DL,VL))



        ptsmat = repeat(pts,1, size(d)[1])
        

#        println("A2 ", A2.(pts))
        D1 = (A1.(ptsmat)) .*(d)

        #D2 = A2(1.0)*(d*d)
        #D1 = A1(1.0)*(d)
        
        v = diagm(A0.(pts))
        #v = A0.(ptsmat)
        
        
        N = size(d)[1] - 1
        
        x = getCpts(N)
        w = baryW1(N)
        y = getCpts1(N-1)
        P = barymat(y, x, w)

        H0 = v + D1
        
        n2r += size(v)[1] - 1
        n2c += size(v)[1] 
        H[(n1r+1):n2r,(n1c+1):n2c] += P*H0
        S[(n1r+1):n2r,(n1c+1):n2c] += P*I(N+1)

        rhs[(n1r+1):n2r] = P*B.(pts)
        
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
