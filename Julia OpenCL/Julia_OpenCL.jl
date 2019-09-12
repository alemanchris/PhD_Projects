#=
 This code uses ArrayFire to Optimize Paralell Computing on the GPU
=#
#--------------------------------#
#         House-keeping          #
#--------------------------------#

using Distributions
using ArrayFire
using Distributed
using Profile
#--------------------------------#
#         Initialization         #
#--------------------------------#
#set JULIA_NUM_THREADS = 1
# Go to platformio-ide-terminal
# cd C:\Users\Aleman\AppData\Local\Julia-1.0.2\bin set JULIA_NUM_THREADS = 1
# Check number of threads
Threads.nthreads()
# Change number of threads

# Grid for x
nx  = 1500
xmin = 0.1
xmax = 4.0

# Grid for e: parameters for Tauchen
ne  = 50 #15
ssigma_eps      = 0.02058
llambda_eps     = 0.99
m               = 1.5

# Utility function
ssigma   = 2
bbeta    = 0.97
T        = 10

# Prices
r  = 0.07
w  = 5

# Initialize the grid for X
xgrid = zeros(1,nx)

# Initialize the grid for E and the transition probability matrix
egrid = zeros(1,ne)
P     = zeros(ne, ne)

# Initialize value function V

V_gpu_b        = zeros(T, nx, ne)
V_tomorrow_gpu_b = zeros(nx, ne)
#--------------------------------#
#         Grid creation          #
#--------------------------------#

# Grid for capital (x)
size_a = nx
xstep = (xmax - xmin) /(size_a - 1)

for i = 1:nx
  xgrid[i] = xmin + (i-1)*xstep
end

# Grid for productivity (e) with Tauchen (1986)
size_e = ne
ssigma_y = sqrt((ssigma_eps^2) / (1 - (llambda_eps^2)))
estep = 2*ssigma_y*m / (size_e-1)
for i = 1:ne
  egrid[i] = (-m*sqrt((ssigma_eps^2) / (1 - (llambda_eps^2))) + (i-1)*estep)
end

# Transition probability matrix (P) Tauchen (1986)
mm = egrid[2] - egrid[1]
function transition_f(ne,egrid,llambda_eps,mm,ssigma_eps)
    for j = 1:ne
          for k = 1:ne
                if (k == 1)
                  P[j, k] = cdf(Normal(), (egrid[k] - llambda_eps*egrid[j] + (mm/2))/ssigma_eps)
                elseif (k == ne)
                  P[j, k] = 1 - cdf(Normal(), (egrid[k] - llambda_eps*egrid[j] - (mm/2))/ssigma_eps)
                else
                  P[j, k] = cdf(Normal(), (egrid[k] - llambda_eps*egrid[j] + (mm/2))/ssigma_eps) - cdf(Normal(), (egrid[k] - llambda_eps*egrid[j] - (mm/2))/ssigma_eps)
                end
          end
    end
    return P
end
P = transition_f(ne,egrid,llambda_eps,mm,ssigma_eps)
# Exponential of the grid e
for i = 1:ne
  egrid[i] = exp(egrid[i])
end
# Vectorize
P_a = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
xgrid_a  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
xgridp_a = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
egrid_a  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])
cons_a = (1+r).*xgrid_a .+  egrid_a.*w .- xgridp_a
neg_ind_a = convert(Array{Bool},cons_a.<=zeros(nx,ne,nx)) # convert(Array{Int}, a)
ut_a = (cons_a.^(1-ssigma))./(1-ssigma)

struct modelState_zx
  ne::Int64
  nx::Int64
  T::Int64
  age::Int64
  P::Array{Float64,2}
  ut_a::Array{Float64,3}#Vector{Float64}
  neg_ind_a::Array{Bool,3}#Vector{Float64}
  ssigma::Float64
  bbeta::Float64
  V::Array{Float64,2}
  w::Float64
  r::Float64
end

function solve_age_gpu_c(currentState::modelState_zx)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P       = currentState.P
    ut_g    = currentState.ut_a
    neg_ind_g   = currentState.neg_ind_a
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    #P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    #V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])
    VV = zeros(nx,ne)
    for ie = 1:ne
        # Calculate expected utility
        #expected = zeros(nx,ne,nx)
        expected = AFArray(zeros(nx,1))
        if (age < T)
            #expected = repeat(sum(P[ie,:].*V,dims=2),outer = [1,1,nx])
            expected = AFArray(reshape(sum(V.*reshape(P[ie,:],(1,ne)),dims=2),(nx,1)))
        end

        # Manipulate Matrices Budget constraint
        #xgrid_g  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
        #xgridp_g = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
        #egrid_g  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])

        # Compute consumption
        #cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
        # Compute Utility
        utility = AFArray(ut_g[:,ie,:]).+ bbeta.*expected
        # Eliminate negatives
        #neg_ind = cons.<=zeros(nx,ne,nx)
        #utility = utility.*(.!reshape(neg_ind_g[:,ie,:],(nx,nx))).+reshape(neg_ind_g[:,ie,:],(nx,nx)).*(-10.0^(5))
        utility = utility.*(AFArray(convert(Array{Bool},.!neg_ind_g[:,ie,:]))).+AFArray(neg_ind_g[:,ie,:].*(-10.0^(5)))
        #VV[:,ie] = Array(reshape(maximum(utility,dims = 2),(nx,1))) #reshape is allowed
        VV[:,ie] = Array(maximum(utility,dims = 2)) #reshape is allowed
    end
    return VV
end

function solve_VFI_gpu_c(T,ne,nx,P,ut_a,neg_ind_a,ssigma,bbeta,w,r,V_tomorrow_gpu_b,V_gpu_b)
    for age = T:-1:1
        currentState = modelState_zx(ne,nx,T,age,P,ut_a,neg_ind_a,ssigma,bbeta,V_tomorrow_gpu_b,w,r)
        tempV_gpu_b = solve_age_gpu_c(currentState)
        V_gpu_b[age,:,:] = tempV_gpu_b
        V_tomorrow_gpu_b = tempV_gpu_b
    end

    return V_gpu_b
end

# GPU
gpu_vector_c = @elapsed V_gpu_solved_c = solve_VFI_gpu_c(T,ne,nx,P,ut_a,neg_ind_a,ssigma,bbeta,w,r,V_tomorrow_gpu_b,V_gpu_b)
