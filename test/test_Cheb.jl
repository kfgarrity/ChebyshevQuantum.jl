using ChebyshevQuantum
using Test

tol_val = 1e-8


@testset "test Cheb interp" begin

    function f(x)
        return exp(-0.2*sin(3*x))
    end

    c = ChebyshevQuantum.Chebyshev.make_Cheb(f , a=2, b = 3)

    @test isapprox(c(2.3), f(2.3), atol=tol_val)
    @test isapprox(c(2.0), f(2.0), atol=tol_val)

    function f2(x)
        return max(sin(5*x), 0.0)
    end

    c2 = ChebyshevQuantum.Chebyshev.make_Cheb(f2, thr=1e-12)

    for x = -1.0:.1:1
        @test isapprox(c2(x), f2(x), atol=tol_val)
    end

    c3 = ChebyshevQuantum.Chebyshev.make_Cheb(f2, thr=1e-12, a=-2, b = 2)

    for x = -2.0:.25:2
        @test isapprox(c3(x), f2(x), atol=tol_val)
    end
    
end

