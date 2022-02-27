

@time begin
@time begin
    function f(x)
        return x
    end


    function y(f, x)
        return f(x)
    end

    @time f(1.0)
    @time y(f, 1.0)
    
#    @time f.([1.0, 2.0])
#    @time y.(f, [1.0, 2.0])
end
end
