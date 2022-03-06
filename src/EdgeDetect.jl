module EdgeDetect

using ForwardDiff
using Plots
using Statistics


function approx_diff(f, x, delta)
    h=1e-6
    return (f(x+delta*h) - f(x-delta*h))/(2*h)
end



function find(f,  N=110; a = -1.0, b = 1.0, deriv=:forwarddiff)

#    println("a")
    
    discont_f0, stddev_f0 = find_edges(f, N, a=a, b = b, deriv=deriv)

#    println("EDGES F0 ", discont_f0)

    begin
        if deriv == :approx
            df = x -> approx_diff(f,x, b-a)
        else
#            println("ForwardDiff")
            df = x -> ForwardDiff.derivative(f, x);
        end
        
        discont_f1, stddev_f1 = find_edges(df, N, a=a, b = b, deriv=deriv)
#        println("EDGES F1 ", discont_f1)
    end
    
    edges = [discont_f0;discont_f1]

    if size(edges)[1] == 0
        return zeros(0,2), Bool[]
    end
        
    #
    #println("EDGES ", edges)
    
    inds = sortperm(edges[:,1])
    edges = edges[inds, :]

    #println("EDGES2 ", edges)

    
    good = [1]
    for i = 2:size(edges)[1]
        if abs(edges[i-1,1] - edges[i,1]) > 1e-3 *(b-a) 
            push!(good,i)
        end
    end
    edges = edges[good,:]
    

    contin = Bool[]
    for i = 1:size(edges)[1]

#        println("contin ", abs(f(edges[i,1]) - f(edges[i,2]) ), " std ", stddev_f0, " r ", abs(f(edges[i,1]) - f(edges[i,2]) ) / stddev_f0 , " r2 ", f(edges[i,1]) / stddev_f0)
        
        if  abs(f(edges[i,1]) - f(edges[i,2]) ) / stddev_f0 < 0.5e-2 && f(edges[i,1]) / stddev_f0 < 1e3
            push!(contin, true)
        else
            push!(contin, false)
        end
    end
    return edges, contin
    
end

function find_edges(f, N = 110; a=-1.0, b = 1.0, deriv=:forwarddiff)



    if deriv == :approx
        df = x -> approx_diff(f,x, b-a)
    else
        df = x -> ForwardDiff.derivative(f, x);
    end


    #    f(0.1)
#    df(0.1)
    grid, grid2, fg, dfg, df_approx, df_err = get_derivs(f,df, N, a=a, b=b)



    
    stddev = std(fg)
    
    local_max = []
    for i in 2:N-1
        if abs(df_err[i])+1e-10 > abs(df_err[i-1]) + 1e-10 && abs(df_err[i]) > abs(df_err[i+1]) + 1e-10
            push!(local_max, i)
        end
    end

#    println("local_max ", local_max, grid[local_max])

    true_edge = zeros(0,2)
    
    for i in local_max
#        println([i, grid[i], grid[i+1], fg[i], fg[i+1], df_err[i] ])

        gridX, grid2X, fgX, dfgX, df_approxX, df_errX = get_derivs(f,df, 20, a=grid[i], b=grid[i+1])
        ind =argmax(abs.(df_errX))
        new_max_err = abs.(df_errX[ind])
        old_max_err = abs(df_err[i])

        if new_max_err >= old_max_err *0.1 || new_max_err > 1e2 #true discont doesn't shrink (much) as we get closer

            for j = 1:7
                gridX, grid2X, fgX, dfgX, df_approxX, df_errX = get_derivs(f,df, 20, a=gridX[ind], b=gridX[ind+1])
                ind =argmax(abs.(df_errX))
                new_max_err = abs.(df_errX[ind])
#                println("iter $j ", [gridX[ind], gridX[ind+1]], " " , gridX[ind+1] - gridX[ind], " " , [new_max_err, old_max_err])
                old_max_err = new_max_err
            end

            if new_max_err >= old_max_err *0.1 && new_max_err > 1e3
                true_edge = vcat(true_edge, [gridX[ind]  gridX[ind+1]])
#                println("TRUE dis ", [gridX[ind], gridX[ind+1]])
#                println()
            end
            
        end

    end

#    println("TRUE EDGE ", true_edge)
    
#    plot(grid2[2:end-1], df_err[2:end-1], label="df err ")
#    plot!(grid2, dfg, label="FD", linewidth=3)
#    plot!(grid2, df_approx, label="df appprox", linewidth=2)
#    plot!(grid, fg, label="f")

    
    return true_edge, stddev
    
end

function get_derivs(f, df, N; a=-1, b = 1)

    
    begin

        aa = a + (b-a)/N * 0.001
        bb = b - (b-a)/N * 0.0015

        h = (bb - aa)/(N) * (1.0 - 10^-12)
        grid = (aa):h:(bb)
        grid2 = (aa:h:bb)[1:end-1] .+ h/2
    end


    
    #    N = length(grid)-2
    
    #println("grid ")

    fg = f.(grid)
    dfg = df.(grid2)
    
    #fg = zeros(length(grid))
    #for i = 1:length(grid)
    #    fg[i] = f(grid[i])
    #end
    #dfg = zeros(length(grid2))    
    #for i = 1:length(grid2)
    #    dfg[i] = df(grid2[i])
    #end
    
#    begin
#        fg = zeros(size(grid))
#        dfg = zeros(size(grid2))
#        for i = 1:length(grid)
#        fg[i] = f(grid[i])
#        end
#        for i = 1:length(grid2)
#            dfg[i] = df(grid2[i])
#        end
#    end    

#    println("size fg ", size(fg))
#    println("N $N")
    
    begin
        df_approx = zeros(N)
        for i in 1:(N)
            
            df_approx[i] = (fg[i+1] - fg[i]) / h
            
            
        end
        df_err = df_approx - dfg
    end

    return grid, grid2, fg, dfg, df_approx, df_err


end    

end #end module
