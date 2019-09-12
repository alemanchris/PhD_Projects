using ArrayFire

maxIterations = 500
gridSize = 2048
xlim = [-0.748766713922161, -0.748766707771757]
ylim = [ 0.123640844894862,  0.123640851045266]

x = range( xlim[1], xlim[2], length = gridSize )
y = range( ylim[1], ylim[2], length = gridSize )

xGrid = [i for i in x, j in y]  # puts x[1] in i,j very funny repmat
yGrid = [j for i in x, j in y]

z0 = xGrid + im*yGrid #im declares the imaginary part (i)

function mandelbrotCPU(z0, maxIterations)
    z = copy(z0)
    count_a = ones( size(z) )

    for n in 1:maxIterations
        z .= z.*z .+ z0
        count_a .+= abs.( z ).<=2
    end
    count_a = log.( count_a )
end

function mandelbrotGPU(z0, maxIterations)
    z = z0
    count_a = ones(AFArray{Float32}, size(z) )

    for n in 1:maxIterations
        z = z .* z .+ z0
        count_a = count_a + (abs(z)<= 2) #mixing AF arrays with non AF Arrays? NO!,
        # cuz he is tricky enough to convert the "z0" input to AFArray before it goes in
    end
    sync(log( count_a )) # Sync to wait until
    # Why doesnt he replace count with the log(count)?, ANSWER: JUST CUZ HE IS DAMM LAZY DAMMIT
end

# warmup
count_a = mandelbrotCPU(z0, 1)

cpu_time = @elapsed count_a = mandelbrotCPU(z0, maxIterations)

count_a .-= minimum(count_a)
count_a ./= maximum(count_a)
img = AFArray(Array{Float32}(count_a))

ArrayFire.image(img) #Rendering an Immage using ArrayFire, I wont use this

count_a = mandelbrotGPU(AFArray(z0), 1)
#count = mandelbrotGPU(AFArray(z0), maxIterations) # This is redundant
gpu_time = @elapsed count_a = mandelbrotGPU(AFArray(z0), maxIterations)

ArrayFire.figure(2)
count_a -= min_all(count_a)[1]
count_a /= max_all(count_a)[1]
img = AFArray{Float32}(count_a) # I could also do AFArray(count)

ArrayFire.image(img)


@show cpu_time, gpu_time, cpu_time/gpu_time
