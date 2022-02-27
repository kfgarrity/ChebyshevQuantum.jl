

@time begin
@time begin
#    function fall()
    function f(x)
        return x
    end


    function y(h, x)
        return h(x)
    end

    println("f")
    @time f(1.0)
    @time f(1.0)
    println("y")
    @time y(f, 1.0)
    @time y(f, 1.0)

    println("f broad")
    @time f.([1.0, 2.0])
    @time f.([1.0, 2.0])
    println("y broad")
    @time y.(f, [1.0, 2.0])
    @time y.(f, [1.0, 2.0])

    function g(x)
        return 2*x
    end

    println("g g broad")
    @time g(1.0)
    @time g.([1.0, 2.0])
    
    
    println("y g broad")
    @time y.(g, [1.0, 2.0])
    @time y.(g, [1.0, 2.0])
    println()
    println("done")
#    end
#
#    @time fall()
#    println("done fall1")
#    @time fall()
#    println("done fall2")
        
    end
    
end
