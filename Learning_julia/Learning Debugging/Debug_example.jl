function enum_mio(n)
    out = zeros(1,n)
    for i =1:n
        out[i]=i
    end
    return out
end

# Tutorial on how to debug.
#=
Juno.@run enum_mio(3) Enter and the continue
Juno.@enter enum_mio(3)
This opens a new GUY and the BAR

=#
