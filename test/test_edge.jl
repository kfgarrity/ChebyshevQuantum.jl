using ChebyshevQuantum
using Test

tol_val = 1e-11


@testset "find edge function " begin

    #step function
    function f(x)
        if x < 0.2
            return 0.0
        else
            return 1.0
        end
    end

    c = ChebyshevQuantum.EdgeDetect.find_edges(f)
    
    @test isapprox(c[1,1], 0.2,  atol=1e-4)
    @test isapprox(c[1,2], 0.2,  atol=1e-4)


    
end


@testset "find edge f, f' " begin

    #two edges
    function f(x)
        if x < 0.3
            return 0.1 * abs(x-0.1) + x^3
        else
            return 1.0 * abs(x-0.1) + x^3
        end
    end

    c = ChebyshevQuantum.EdgeDetect.find(f)
    
    @test isapprox(c[1,1], 0.1,  atol=1e-4)
    @test isapprox(c[1,2], 0.1,  atol=1e-4)
    
    @test isapprox(c[2,1], 0.3,  atol=1e-4)
    @test isapprox(c[2,2], 0.3,  atol=1e-4)


    
end


@testset "find edge singularity " begin

    #singularity
    function f(x)
        return 2.0 * sqrt(abs(x-0.1))
    end

    c = ChebyshevQuantum.EdgeDetect.find(f)
    
    @test isapprox(c[1,1], 0.1,  atol=1e-4)
    @test isapprox(c[1,2], 0.1,  atol=1e-4)

    #singularity
    function f2(x)
        return -max(1.0 / (x-0.1), -1e10)
    end

    c2 = ChebyshevQuantum.EdgeDetect.find(f2)
    
    @test isapprox(c2[1,1], 0.1,  atol=1e-4)
    @test isapprox(c2[1,2], 0.1,  atol=1e-4)
    
    
end

@testset "don't find  " begin

    #no edges
    function f(x)
        return sin(30*x)
    end

    c = ChebyshevQuantum.EdgeDetect.find(f)
    
    @test size(c)[1] == 0

    #no edges
    function f2(x)
        return exp(-1/x^2)
    end

    c2 = ChebyshevQuantum.EdgeDetect.find(f2)
    
    @test size(c2)[1] == 0

    #edge at end, don't find
    function f3(x)
        return max(-1/x, -1e10)
    end

    c3 = ChebyshevQuantum.EdgeDetect.find(f3, a=0, b =1)
    
    @test size(c3)[1] == 0
    
    
end

