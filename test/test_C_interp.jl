using ChebyshevQuantum
using Test

tol_val = 1e-11


@testset "test C interp" begin

    function f(x)
        return exp(-0.2*sin(3*x))
    end

    c = ChebyshevQuantum.Chebyshev.make_Cheb(f , a=2, b = 3)

    @test isapprox(c(2.3), f(2.3), atol=tol_val)

    function f2(x)
        if x < 0.1
            return exp(-0.2*sin(3*x))
        else
            return 1+x^2
        end
    end

    c = ChebyshevQuantum.Chebyshev.make_Cheb(f2 , a=-1, b = 1)

    x=0.0
    @test isapprox(c(x), f2(x), atol=tol_val)
    x=0.2
    @test isapprox(c(x), f2(x), atol=tol_val)
    


end

