# Paralell Computing in Julia
The codes in the following repository contain examples of paralell computing in Julia.

The following codes solve the Life Cycle Model example presented on: A Practical Guide to Parallelization in Economics by Villaverde & Zarruk (2018)
- Julia_Par.jl uses the Julia Distributed library to paralellize over the CPU (this is taken directly from Villaverde & Zarruk(2018)).
- Julia_OpenCL.jl uses ArrayFire for Julia to parallelize over the GPU.

UPDATE COMMING: Now I have NVIDIA GPU, stay tuned for the code written for Julia using CUDA.

The code is written for Julia 1.0.2 (Sorry for not keeping up with the latest version)

Make sure to have ArrayFire loaded and working on your operating system.

This is the link to Villaverde & Zarruk(2018) repository

https://github.com/davidzarruk/Parallel_Computing
