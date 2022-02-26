using ChebyshevQuantum
using LinearAlgebra



function go(N)

    D = ChebyshevQuantum.Interp.getD(N);
    pts = ChebyshevQuantum.Interp.getCpts(N);
    #v = eigvals( -(D*D )[2:N,2:N]);
    v, _ = eigs( -(D*D )[2:N,2:N], nev=5, which=:SR, tol=1e-4, maxiter=30000);
    v / (pi^2/4)
end

function go_dense(N)

    D = ChebyshevQuantum.Interp.getD(N);
    pts = ChebyshevQuantum.Interp.getCpts(N);
    v = eigvals( -(D*D )[2:N,2:N]);
    #v, _ = eigs( -(D*D )[2:N,2:N], nev=10, which=:SR);
    v / (pi^2/4)
end

