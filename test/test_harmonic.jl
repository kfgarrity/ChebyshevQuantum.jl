using ChebyshevQuantum
using Test

tol_val = 1e-8


@testset "test harmonic ω=50" begin


    ω = 50.0
    
    function V(x)
        return 0.5 * ω^2 * x^2
    end
    
    vals, vects = ChebyshevQuantum.SolveEig.solve(V=V, N=60, dosplit=false)
    #H = ChebyshevQuantum.SolveEig.setup_ham(V, 60, a=-1.0, b= 1.0);
    #vals, vects =  ChebyshevQuantum.SolveEig.solve(H);

    @test isapprox(vals[1],   (0.5+0) * ω)
    @test isapprox(vals[2],   (0.5+1) * ω)
    @test isapprox(vals[3],   (0.5+2) * ω)
    @test isapprox(vals[4],   (0.5+3) * ω)
    @test isapprox(vals[5],   (0.5+4) * ω)
    
end


@testset "test harmonic ω=1, " begin

    ω = 1.0
    
    function V(x)
        return 0.5 * ω^2 * x^2
    end
    

    vals, vects = ChebyshevQuantum.SolveEig.solve(V=V, N=60, a=-10.0, b= 10.0, dosplit=false)

    #    H = ChebyshevQuantum.SolveEig.setup_ham(V, 60, a=-10.0, b= 10.0);
    #    vals, vects =  ChebyshevQuantum.SolveEig.solve(H);

    @test isapprox(vals[1],   (0.5+0) * ω)
    @test isapprox(vals[2],   (0.5+1) * ω)
    @test isapprox(vals[3],   (0.5+2) * ω)
    @test isapprox(vals[4],   (0.5+3) * ω)
    @test isapprox(vals[5],   (0.5+4) * ω)
    
end
