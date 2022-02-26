using ChebyshevQuantum
using LinearAlgebra



function go(N ; Rmax = 20, L=0)

    @time begin 
        pts = ChebyshevQuantum.Interp.getCpts(N);
        D = ChebyshevQuantum.Interp.getD(pts);
    end

    @time begin
        r = (1 .+ pts) * Rmax / 2
        H1 = -0.5*(D*D )[2:N,2:N] * 4/Rmax^2
        H2 = -diagm(1 ./ r[2:N])
        H3 = L*(L+1)*diagm(1 ./ r[2:N].^2) / 2.0
        H = H1+H2+H3
    end
    
    @time v = eigvals( H )

    m = min(10, length(v))
    return v[1:m] * 2 .* ( (1:m) .+ L).^2
end


function go_exp(N ; Rmax = 20, L=0, a = 1.0)

    @time begin 
        pts = ChebyshevQuantum.Interp.getCpts(N);
        D = ChebyshevQuantum.Interp.getD(pts);
    end

    @time begin
        r = (1 .+ pts) * Rmax / 2
        H1 = -0.5*(D*D )[2:N,2:N] * 4/Rmax^2 

        #H1a = zeros(N-1,N-1)
        H1a = (D )[2:N,2:N] * 2/Rmax * a

        #        H1b = -0.5*I(N-1) * a^2
        #H1b = zeros(N-1, N-1)
        
        H2 = -diagm(1 ./ r[2:N])
        H3 = L*(L+1)*diagm(1 ./ r[2:N].^2) / 2.0
        H = H1+H1a+H2+H3
    end
    
    @time v = eigvals( H )
    v = v .- 0.5 * a^2
    m = min(10, length(v))
    return v[1:m] * 2 .* ( (1:m) .+ L).^2
end

#=function go_exp2(N ; Rmax = 20, L=0, a = 1.0)

    @time begin 
        pts = ChebyshevQuantum.Interp.getCpts(N);
        D = ChebyshevQuantum.Interp.getD(pts);
    end

    @time begin
        r = (1 .+ pts) * Rmax / 2
        H1 = -0.5*(D*D )[2:N,2:N] * 4/Rmax^2 

        #H1a = zeros(N-1,N-1)
        H1a = (D )[2:N,2:N] * 2/Rmax * a

        H1b = -0.5*I(N-1) * a^2
        #H1b = zeros(N-1, N-1)
        
        H2 = -diagm(1 ./ r[2:N])
        H3 = L*(L+1)*diagm(1 ./ r[2:N].^2) / 2.0
        H = H1+H1a+H1b+H2+H3
    end
    
    @time v = eigvals( H )
#    v = v .- 0.5 * a^2
    m = min(10, length(v))
    return v[1:m] * 2 .* ( (1:m) .+ L).^2
end
=#


