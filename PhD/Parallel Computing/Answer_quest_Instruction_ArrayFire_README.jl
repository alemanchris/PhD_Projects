using Random, ArrayFire
#=
 QUESTION 1
"." opperation on CPU vs GPU:
on CPU looks like log.(A) is the log of each element on A; putting log(A) wont work.
on GPU looks like log(A) is allowed for the log of each element
=#
size_a = 3000
A = rand(size_a,size_a)
# Will log(A) work given that it is in the CPU?
#A_1 = log(A) # Se cuelga (THIS SHIT DOESNT WORK: PROVED)
# What is the difference between log.(A)?
A_2 = log.(A)
A_gpu = AFArray(A)
A_gpu2 = log(A_gpu)
A_gpu2 == A_2

A_2[floor(Int,1.0),1] # This is my solution
A_2[convert(Int8,floor(1.0)),1]

# This might answer Question 3
what_type_A = A_gpu2 +1
what_type_B = A_gpu2 + A_2 # I cannot mix AFArrays
what_type_C = A_gpu2 + A_gpu2
what_type_D = A_gpu2 + 1.2
# By mixing types I get an AFfire type
# Proved: log(A) is allowed for GPU arrays but not for CPU arrays. For CPU Arrays use, log.(A)

#=
QUESTION 2
Is the output of a Function that outputs am AFArray also an AF Array?
Apparently yes: if u try to use a non wrapped function on the output it should fail.
=#
function my_test(A_imput)
    B_gpu = AFArray(A_imput)+AFArray(ones(size(A_imput)))
    return B_gpu
end
B_gpu = my_test(A_2)
#=
YES THE ANSWER IS AN AFARRAY
QUESTION 3:
What happens if I mix AFArray computation with a Non AFArray? Will it still be performed?
If YES, where in the CPU or the GPU.
# This might answer Question 3
what_type_A = A_gpu2 +1    # YES if integers
what_type_B = A_gpu2 + A_2 # I cannot mix AFArrays
what_type_C = A_gpu2 + A_gpu2
=#
QUESTION 4:
Can I perform soms stuff that are not in the wrap.jl of the ArrayFire?
# Do some permute
# YES I can do some reshape and the result is a AFArray.
A_3D = reshape(A_2,(size_a,size_a,1))
A_gpu_3D = reshape(A_gpu,(size_a,size_a,1))

# Repeat
A_3D_big = repeat(A_3D,outer = [1,1,2])

# Out of memory, for GPU on a 3d with 3000
# Error if you rin any of these two
# Repeat is not fucked only it takes damm a lot of time
A_gpu_3D_big = repeat(A_gpu_3D,outer = [1,1,2]) #repeat is fucked with AFArray
A_gpu_3D_big = repeat(A_gpu_3D,(1,1,2)) #repeat is fucked with AFArray
# Permute dimensions
# Not cut off for AFArray
A_gpu_permute = permutedims(A_gpu_3D,[2,3,1])
A_permute = AFArray(permutedims(A_3D,[2,3,1]))
A_permute2 = permutedims(A_3D,[2,3,1])
neg_ind = A_permute2.<=(ones(size_a,1,size_a).*(-0.5))
neg_ind_gpu = A_permute.<=AFArray((ones(size_a,1,size_a).*(-0.5)))

QUESTION 5:
Can I mix non AFArrays and AFArrays in the inputs of a function?
    function test_sum(A,B)

        return C
    end
