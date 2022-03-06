module ChebyshevQuantum

using Plots
using ForwardDiff
using Polynomials
using LinearAlgebra

greet() = print("Hello World!")

include("Interp.jl")
include("EdgeDetect.jl")
include("Cheb.jl")
include("Op.jl")
include("SolveEig.jl")
include("SolveDiff.jl")

end # module
