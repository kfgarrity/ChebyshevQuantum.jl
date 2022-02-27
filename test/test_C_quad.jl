using ChebyshevQuantum
using Test


@testset "test C quad and derivative" begin

    function f(x)
        return exp(1.2*x)
    end

    c = ChebyshevQuantum.Chebyshev.make_Cheb(f , a=2, b = 3)

    exact = (exp(1.2*3) - exp(1.2*2))/1.2
    int = ChebyshevQuantum.Chebyshev.integrateCC(c)
    int2 = sum(c)

    @test isapprox(exact, int, atol=1e-9)
    @test isapprox(exact, int2, atol=1e-9)

    x=2.2
    exact = 1.2 * exp(x*1.2)
    
    cprime = ChebyshevQuantum.Interp.derivative(c)
    cprime2 = c'

    @test isapprox(exact, cprime(x), atol=1e-9)
    @test isapprox(exact, cprime2(x), atol=1e-9)
    

end
