using ChebyshevQuantum
using Test

tol_val = 1e-11


@testset "test interp" begin

    function f(x)
        return exp(-0.2*sin(3*x))
    end

    c = ChebyshevQuantum.Interp.make_cheb(f , a=2, b = 3)

    @test isapprox(c(2.3), f(2.3), atol=tol_val)
    @test isapprox(c(2.0), f(2.0), atol=tol_val)
    @test isapprox(c(2.712), f(2.712), atol=tol_val)

    function f2(x)
        return log(x) + x^2
    end

    c2 = ChebyshevQuantum.Interp.make_cheb(f2 , a=2, b = 3)

    @test isapprox(c2(2.3), f2(2.3), atol=tol_val)
    @test isapprox(c2(2.0), f2(2.0), atol=tol_val)
    @test isapprox(c2(2.712), f2(2.712), atol=tol_val)
    


end

@testset "test functions var N" begin

    function f1(x)
        return exp(-0.2*sin(3*x))
    end
    function f2(x)
        return exp(-0.3*cos(3*x)) + 20
    end

    c1 = ChebyshevQuantum.Interp.make_cheb(f1 , a=2, b = 3)
    c2 = ChebyshevQuantum.Interp.make_cheb(f2 , a=2, b = 3)

    x=2.712
    @test isapprox( (c1+c2)(x), f1(x)+f2(x) , atol=tol_val)
    @test isapprox( (c1-c2)(x), f1(x)-f2(x) , atol=tol_val)
    @test isapprox( (c1^2)(x), f1(x)^2 , atol=tol_val)
    @test isapprox( (c1/c2)(x), f1(x)/f2(x) , atol=tol_val)

    @test isapprox( (-c1)(x), -f1(x) , atol=tol_val)
    @test isapprox( (2*c1)(x), 2*f1(x) , atol=tol_val)
    @test isapprox( (c1*2)(x), 2*f1(x) , atol=tol_val)
    @test isapprox( (2+c1)(x), 2+f1(x) , atol=tol_val)
    @test isapprox( (c1+2)(x), f1(x)+2 , atol=tol_val)
    
    @test isapprox( (2 / c2)(x), 2 / f2(x) , atol=tol_val)
    @test isapprox( (2 - c2)(x), 2 - f2(x) , atol=tol_val)

    
end



@testset "test functions fixed N" begin

    function f1(x)
        return exp(-0.2*sin(3*x))
    end
    function f2(x)
        return exp(-0.3*cos(3.5*x)) + 20
    end

    c1 = ChebyshevQuantum.Interp.make_cheb(f1 , a=2, b = 3, N=20, thr=-1)
    c2 = ChebyshevQuantum.Interp.make_cheb(f2 , a=2, b = 3, N=20, thr=-1)

    x=2.712
    @test isapprox( (c1+c2)(x), f1(x)+f2(x) , atol=tol_val)
    @test isapprox( (c1-c2)(x), f1(x)-f2(x) , atol=tol_val)
    @test isapprox( (c1^2)(x), f1(x)^2 , atol=tol_val)
    @test isapprox( (c1/c2)(x), f1(x)/f2(x) , atol=tol_val)

    @test isapprox( (-c1)(x), -f1(x) , atol=tol_val)
    @test isapprox( (2*c1)(x), 2*f1(x) , atol=tol_val)
    @test isapprox( (c1*2)(x), 2*f1(x) , atol=tol_val)
    @test isapprox( (2+c1)(x), 2+f1(x) , atol=tol_val)
    @test isapprox( (c1+2)(x), f1(x)+2 , atol=tol_val)


    @test isapprox( (2 / c2)(x), 2 / f2(x) , atol=tol_val)
    @test isapprox( (2 - c2)(x), 2 - f2(x) , atol=tol_val)
    

    
end

@testset "function returns complex number " begin

    function f(x)
        return exp(-1.0im*x)
    end

    c = ChebyshevQuantum.Interp.make_cheb(f , a=-pi, b = pi)

    x= pi/8
    @test isapprox( real(c(x)), real(f(x)) , atol=tol_val)
    @test isapprox( imag(c(x)), imag(f(x)) , atol=tol_val)

    x= 0.1234
    @test isapprox( real(c(x)), real(f(x)) , atol=tol_val)
    @test isapprox( imag(c(x)), imag(f(x)) , atol=tol_val)

    x= 0.5678
    @test isapprox( real(c)(x), real(f(x)) , atol=tol_val)
    @test isapprox( imag(c)(x), imag(f(x)) , atol=tol_val)
    
end
