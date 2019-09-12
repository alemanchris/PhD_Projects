using ArrayFire
using Test
a = rand(10,10)
ad = AFArray(a)
bd = ad+1
function der_sum!(bd)
    for i=1:5000000
        bd = bd+1
    end
    return(bd)
end
ab = der_sum!(bd)
b = Array(bd)
###
function parallel_add!(y, x)
    Threads.@threads for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end
N = 2^20
x = fill(1.0f0, N)  # a vector filled with 1.0 (Float32)
y = fill(2.0f0, N)
fill!(y, 2)
parallel_add!(y, x)
@test all(y .== 3.0f0)
