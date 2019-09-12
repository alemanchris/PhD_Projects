#using ArrayFire
function enum_mio(n)
    out = zeros(1,n)
#    fire_1 = AFArray(out)
    for i =1:n

        out[i]=i
#        fire_1[i] = i
    end
    return out
end

# Tutorial on how to debug.
#ON THE REPL
#=
Juno.@run enum_mio(3) Enter and the continue
Juno.@enter enum_mio(3)
This opens a new GUY and the BAR

=#
# Does the debugger work with GPU?
