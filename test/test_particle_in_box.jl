using ChebyshevQuantum
using Test

tol_val = 1e-8


@testset "test particle in box " begin

    H = ChebyshevQuantum.SolveEig.setup_ham(x->0.0, 40);
    vals, vects =  ChebyshevQuantum.SolveEig.solve(H)

    @test isapprox(vals[1],   pi^2/8)
    @test isapprox(vals[2], 2^2*pi^2/8)
    @test isapprox(vals[3], 3^2*pi^2/8)
    @test isapprox(vals[4], 4^2*pi^2/8)
    @test isapprox(vals[5], 5^2*pi^2/8)
    
end

@testset "test particle in bigger box " begin

    H = ChebyshevQuantum.SolveEig.setup_ham(x->0.0, 40, a = 0, b = 10);
    vals, vects =  ChebyshevQuantum.SolveEig.solve(H)

    f = (2/10)^2
    
    @test isapprox(vals[1],   pi^2/8 * f)
    @test isapprox(vals[2], 2^2*pi^2/8 * f)
    @test isapprox(vals[3], 3^2*pi^2/8 * f)
    @test isapprox(vals[4], 4^2*pi^2/8 * f)
    @test isapprox(vals[5], 5^2*pi^2/8 * f)
    
end

tol_val = 1e-10

@testset "test particle in box N=100 " begin

    H = ChebyshevQuantum.SolveEig.setup_ham(x->0.0, 100);
    vals, vects =  ChebyshevQuantum.SolveEig.solve(H)

    @test isapprox(vals[1],   pi^2/8)
    @test isapprox(vals[2], 2^2*pi^2/8)
    @test isapprox(vals[3], 3^2*pi^2/8)
    @test isapprox(vals[4], 4^2*pi^2/8)
    @test isapprox(vals[5], 5^2*pi^2/8)
    
end
