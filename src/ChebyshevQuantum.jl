module ChebyshevQuantum
using Plots

greet() = print("Hello World!")

include("Interp.jl")
include("EdgeDetect.jl")
include("Cheb.jl")

end # module
