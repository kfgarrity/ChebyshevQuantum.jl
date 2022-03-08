using ChebyshevQuantum
using Test

tol_val = 1e-10


@testset "test particle in box " begin

    H = ChebyshevQuantum.SolveEig.setup_ham(x->0.0, 40);
    #    vals, vects =  ChebyshevQuantum.SolveEig.solve(H)
    vals, vects = eigen( H[2:end-1,2:end-1]);
    

    @test isapprox(vals[1],   pi^2/8, atol=tol_val)
    @test isapprox(vals[2], 2^2*pi^2/8, atol=tol_val)
    @test isapprox(vals[3], 3^2*pi^2/8, atol=tol_val)
    @test isapprox(vals[4], 4^2*pi^2/8, atol=tol_val)
    @test isapprox(vals[5], 5^2*pi^2/8, atol=tol_val)
    
end

@testset "test particle in bigger box " begin

    H = ChebyshevQuantum.SolveEig.setup_ham(x->0.0, 40, a = 0, b = 10);
#    vals, vects =  ChebyshevQuantum.SolveEig.solve(H)
    vals, vects = eigen( H[2:end-1,2:end-1]);
    

    f = (2/10)^2
    
    @test isapprox(vals[1],   pi^2/8 * f, atol=tol_val)
    @test isapprox(vals[2], 2^2*pi^2/8 * f, atol=tol_val)
    @test isapprox(vals[3], 3^2*pi^2/8 * f, atol=tol_val)
    @test isapprox(vals[4], 4^2*pi^2/8 * f, atol=tol_val)
    @test isapprox(vals[5], 5^2*pi^2/8 * f, atol=tol_val)
    
end

tol_val = 1e-10

@testset "test particle in box N=100 " begin

    H = ChebyshevQuantum.SolveEig.setup_ham(x->0.0, 100);
#    vals, vects =  ChebyshevQuantum.SolveEig.solve(H)
    vals, vects = eigen( H[2:end-1,2:end-1]);

    @test isapprox(vals[1],   pi^2/8, atol=tol_val)
    @test isapprox(vals[2], 2^2*pi^2/8, atol=tol_val)
    @test isapprox(vals[3], 3^2*pi^2/8, atol=tol_val)
    @test isapprox(vals[4], 4^2*pi^2/8, atol=tol_val)
    @test isapprox(vals[5], 5^2*pi^2/8, atol=tol_val)
    
end

tol_val = 1e-9

@testset "test particle in finite box " begin

    #test edge detection
    
    function VV(x)
        if abs(x) >= 1.0
            return 10.0
        else
            return 0.0
        end
    end


#    H, S = ChebyshevQuantum.SolveEig.total_ham_piecewise(VV, 30, [-2.0, -1.0, 1.0, 2.0], d2=-0.5)
    #    H, S, lhs = ChebyshevQuantum.Op.setup(A2=-0.5, A0=VV,a=-2.0,b=2.0, N=30)

    vals, vects = ChebyshevQuantum.SolveEig.solve(V=VV, a=-2.0, b = 2.0)

    #    vals  =  eigvals(H, S)
#    vals = real.(vals)
    
    @test isapprox(vals[1],  0.819848145939759, atol=tol_val) #reference values from chebfun.m (see below)
    @test isapprox(vals[2],  3.222122379353262, atol=tol_val)

end

"""
matlab code using chebfun
tic
x = chebfun('x', [-2 2]);
V = 10.0*(abs(x)>=1);
L = chebop(-2,2);
L.op = @(x,u) -0.5*diff(u,2) + V*u;
L.bc = 0;
neigs = 12;
[EV,D] = eigs(L,neigs);
disp(diag(D))
toc
"""


nothing
