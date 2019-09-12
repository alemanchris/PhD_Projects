using ArrayFire
using BenchmarkTools
#ad = AFArray(a)
N = 2^25
JULIA_NUM_THREADS=4
#test_vaue = AFArray(N)
#function mfu_gpu()
#end
function parallel_add!(y, x)
    Threads.@threads for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end
function sequential_add!(y, x)
    for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end
    return nothing
end
function GPU_add!(y, x)
    y = y+ x
    return nothing
end
# Serial
a_1 = fill(1.0f0, N)  # a vector filled with 1.0 (Float32)
b_1 = fill(2.0f0, N)
@btime sequential_add!($a_1, $b_1)
@elapsed b_1 = b_1+a_1
# CPU
c_1 = fill(1.0f0, N)  # a vector filled with 1.0 (Float32)
d_1 = fill(2.0f0, N)
@btime parallel_add!($c_1, $d_1)
# GPU
x = AFArray(fill(1.0f0, N))  # a vector filled with 1.0 (Float32)
y = AFArray(fill(2.0f0, N))
@elapsed y = y+ x
