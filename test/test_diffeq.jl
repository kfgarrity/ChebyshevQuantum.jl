using ChebyshevQuantum
using Test
using SpecialFunctions

tol_var=1e-8

@testset "bessel eq test" begin


    a=0.0;
    b= 5.520078110588564; #second zero of besselj0

    #solve bessels equation between 0 and 2nd zero of besselj0 with Direchlet b.c.'s.
    t = ChebyshevQuantum.DiffEq.solve(bc1 = [:a,0, 1.0], bc2 = [:b,0,0.0], A0=x->x^2, A1= x->x, A2 = x->(x^2) ,B = 0.0, N=20,  a=a, b= b);

    x = a:.01:b

    @test maximum(abs.(besselj0.(x) - t.(x))) < tol_var
    @test maximum(abs.(x.^2 .* t''.(x) + x .* t'.(x) + (x.^2 ).*t.(x))) < tol_var #test diffeq directly

                  
    a=1.0;
 #   b= 3.8319091796875; #first zero of besselj1
    b=   3.8317059218883514
    #solve bessels equation α = 1.0  between 1 and zero of besselj1 with Direchlet b.c.'s.
    t = ChebyshevQuantum.DiffEq.solve(bc1 = [:a,0, besselj1(a)], bc2 = [:b,0,besselj1(b)], A0=x->(x^2 - 1.0), A1= x->x, A2 = x->x^2  ,B = 0.0, N=20,  a=a, b= b);

    x = a:.01:b

    @test maximum(abs.(besselj1.(x) - t.(x))) < tol_var
    @test maximum(abs.(x.^2 .* t''.(x) + x .* t'.(x) + (x.^2 .- 1 ).*t.(x))) < tol_var #test diffeq directly


    
end



@testset "heat eq test" begin


    a=0.0;
    b= 1.0; 

    #solve heat equation
    t = ChebyshevQuantum.DiffEq.solve(bc1 = [:a,0,1.0], bc2 = [:b,0,0.5], A2 = 0.5 , N=20,  a=a, b= b);

    x = a:.01:b

    function f(x) #solution is a line
        1.0 - 0.5*x
    end
    
    @test maximum(abs.(f.(x) - t.(x))) < tol_var
    

end

@testset "SHO test" begin


    a=0.0;
    b=2*pi; 

    #solve simple harmonic oscillator equation, neuman b.c.'s
    t = ChebyshevQuantum.DiffEq.solve(bc1 = [:a,0,1.0], bc2 = [:a,1,0.0], A2 = 1.0, A0=1.0 , N=20,  a=a, b= b);

    x = a:.01:b

    function f(x) #solution is a cosine
        cos(x)
    end
    
    @test maximum(abs.(f.(x) - t.(x))) < tol_var
    

end


@testset "Decay Exp test" begin


    a=0.0;
    b=2*pi; 

    #solve decaying exp equation, neuman b.c.'s
    t = ChebyshevQuantum.DiffEq.solve(bc1 = [:a,0,1.0], bc2 = [:a,1,-1.0], A2 = 1.0, A0=-1.0 , N=20,  a=a, b= b);

    x = a:.01:b

    function f(x) #solution is a decaying exp
        exp(-x)
    end
    
    @test maximum(abs.(f.(x) - t.(x))) < tol_var
    

end

@testset "First order Decay Exp" begin


    a=0.0;
    b=10.0; 

    #solve decaying exp equation, first order
    t = ChebyshevQuantum.DiffEq.solve(bc1 = [:a,0,1.0],  A1 = 1.0, A0=1.0 , N=20,  a=a, b= b);

    x = a:.01:b

    function f(x) #solution is a decaying exp
        exp(-x)
    end
    
    @test maximum(abs.(f.(x) - t.(x))) < tol_var
    @test maximum(abs.(  t'.(x) .+ t.(x))) < tol_var #test diffeq directly
    

end


