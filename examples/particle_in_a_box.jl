using ChebyshevQuantum
using LinearAlgebra



function go(N)

    D = ChebyshevQuantum.Interp.getD(N);
    pts = ChebyshevQuantum.Interp.getCpts(N);
    v = eigvals( (D*D )[2:N,2:N]);
    v / (pi^2/4)
end
