module ChebyshevQuantum

using Plots
using ForwardDiff
using Polynomials
using LinearAlgebra

greet() = print("Hello World!")

include("Interp.jl")
include("EdgeDetect.jl")
include("Cheb.jl")
include("SolveEig.jl")

end # module
