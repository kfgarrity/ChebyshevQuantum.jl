using ChebyshevQuantum
using Test


@testset "test interp" begin

    function f(x)
        return exp(-0.2*sin(3*x))
    end

    c = ChebyshevQuantum.Interp.make_cheb(f , a=2, b = 3)

    @test isapprox(c(2.3), f(2.3), atol=1e-9)
    @test isapprox(c(2.0), f(2.0), atol=1e-9)
    @test isapprox(c(2.712), f(2.712), atol=1e-9)


    x=2.712
    @test isapprox(c(x)+1, (c+1)(x), atol=1e-9)
    @test isapprox(c(x)-1, (c-1)(x), atol=1e-9)

    x=2.0
#    println(c(x), " " , f(x))
#    println(c(x)*2)
#    println( (c*2)(x))

#    println()
#    println(c.f)
#    println((c*2).f)

             
    @test isapprox(c(x)*1.2, (c*1.2)(x), atol=1e-9)
    @test isapprox(c(x)/1.2, (c/1.2)(x), atol=1e-9)

end
