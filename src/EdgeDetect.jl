module EdgeDetect

using ForwardDiff
using Plots

function find(f::Function, N=110; a = -1.0, b = 1.0)

    println("discont in f")
    discont_f0 = find_edges(f, N, a=a, b = b)
    println()
    println()
    println("discont in f prime ")    
    df = x -> ForwardDiff.derivative(f, x);
    discont_f1 = find_edges(df, N, a=a, b = b)

    edges = [discont_f0;discont_f1]
    return edges
    
end

function find_edges(f::Function, N = 110; a=-1.0, b = 1.0)

    df = x -> ForwardDiff.derivative(f, x);

    grid, grid2, fg, dfg, df_approx, df_err = get_derivs(f,df, N, a=a, b=b)
    
    local_max = []
    for i in 2:N-1
        if abs(df_err[i])+1e-10 > abs(df_err[i-1]) + 1e-10 && abs(df_err[i]) > abs(df_err[i+1]) + 1e-10
            push!(local_max, i)
        end
    end

#    println("local_max ", local_max)

    true_edge = zeros(0,2)
    
    for i in local_max
#        println([i, grid[i], grid[i+1], fg[i], fg[i+1], df_err[i] ])

        gridX, grid2X, fgX, dfgX, df_approxX, df_errX = get_derivs(f,df, 20, a=grid[i], b=grid[i+1])
        ind =argmax(abs.(df_errX))
        new_max_err = abs.(df_errX[ind])
        old_max_err = abs(df_err[i])
        if new_max_err >= old_max_err *0.5 #true discont doesn't shrink (much) as we get closer

            for j = 1:4
                gridX, grid2X, fgX, dfgX, df_approxX, df_errX = get_derivs(f,df, 20, a=gridX[ind], b=gridX[ind+1])
                ind =argmax(abs.(df_errX))
                new_max_err = abs.(df_errX[ind])
                println("iter $j ", [gridX[ind], gridX[ind+1]], " " , gridX[ind+1] - gridX[ind], " " , [new_max_err, old_max_err])
                old_max_err = new_max_err
            end

            if new_max_err >= old_max_err *0.5 && new_max_err > 1e4
                true_edge = vcat(true_edge, [gridX[ind]  gridX[ind+1]])
                println("TRUE dis ", [gridX[ind], gridX[ind+1]])
#                println()
            end
            
        end

    end

    println("TRUE EDGE ", true_edge)
    
    plot(grid2[2:end-1], df_err[2:end-1], label="df err ")
    plot!(grid2, dfg, label="FD", linewidth=3)
    plot!(grid2, df_approx, label="df appprox", linewidth=2)
    plot!(grid, fg, label="f")

    return true_edge
    
end

function get_derivs(f, df, N; a=-1, b = 1)

    h = (b-a)/N * (1.0 - 10^-12)
    grid = a:h:b
    grid2 = (a:h:b)[1:end-1] .+ h/2

    fg = f.(grid)

    dfg = df.(grid2)
    
    df_approx = zeros(N)
    for i in 1:(N)

        df_approx[i] = (fg[i+1] - fg[i]) / h

        
    end
    df_err = df_approx - dfg
    return grid, grid2, fg, dfg, df_approx, df_err


end    

end #end module
