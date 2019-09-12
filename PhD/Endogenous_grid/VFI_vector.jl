
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
V_cpu          = zeros(T, nx, ne)
V_gpu          = zeros(T, nx, ne)
V_gpu_a        = zeros(T, nx, ne)
V_gpu_b        = zeros(T, nx, ne)
V_tomorrow_cpu = zeros(nx, ne)
V_tomorrow_gpu = zeros(nx, ne)#AFArray(zeros(nx, ne))
V_tomorrow_gpu_a = zeros(nx, ne)
V_tomorrow_gpu_b = zeros(nx, ne)
# Initialize value function as a shared array
tempV_cpu = zeros(ne*nx) # Used for thread guy
#tempV_gpu = zeros(nx,ne) #Not used
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
struct modelState_a
  ne::Int64
  nx::Int64
  T::Int64
  age::Int64
  P::Array{Float64,2}
  xgrid::Array{Float64,2}#Vector{Float64}
  egrid::Array{Float64,2}#Vector{Float64}
  ssigma::Float64
  bbeta::Float64
  V::Array{Float64,2}
  w::Float64
  r::Float64
end

struct modelState_b
  ind::Int64
  ne::Int64
  nx::Int64
  T::Int64
  age::Int64
  P::Array{Float64,2}
  xgrid::Array{Float64,2}#Vector{Float64}
  egrid::Array{Float64,2}#Vector{Float64}
  ssigma::Float64
  bbeta::Float64
  V::Array{Float64,2}
  w::Float64
  r::Float64
end

struct modelState_x
  ne::Int64
  nx::Int64
  T::Int64
  age::Int64
  P_a::Array{Float64,3}
  xgrid_a::Array{Float64,3}#Vector{Float64}
  xgridp_a::Array{Float64,3}#Vector{Float64}
  egrid_a::Array{Float64,3}#Vector{Float64}
  ssigma::Float64
  bbeta::Float64
  V::Array{Float64,2}
  w::Float64
  r::Float64
end

struct modelState_y
  ne::Int64
  nx::Int64
  T::Int64
  age::Int64
  P_b::Array{Float64,3}
  ut_b::Array{Float64,3}#Vector{Float64}
  neg_ind_b::BitArray{3}#Vector{Float64}
  ssigma::Float64
  bbeta::Float64
  V::Array{Float64,2}
  w::Float64
  r::Float64
end

struct modelState_yy
  ne::Int64
  nx::Int64
  T::Int64
  age::Int64
  P_b::Array{Float64,2}
  ut_b::Array{Float64,3}#Vector{Float64}
  neg_ind_b::BitArray{3}#Vector{Float64}
  ssigma::Float64
  bbeta::Float64
  V::Array{Float64,2}
  w::Float64
  r::Float64
end

struct modelState_zz
  ne::Int64
  nx::Int64
  T::Int64
  age::Int64
  P_a::Array{Float64,3}
  ut_a::Array{Float64,3}#Vector{Float64}
  neg_ind_a::Array{Bool,3}#Vector{Float64}
  ssigma::Float64
  bbeta::Float64
  V::Array{Float64,2}
  w::Float64
  r::Float64
end

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

P_c = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
xgrid_c  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
xgridp_c = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
egrid_c  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])
cons_c = (1+r).*xgrid_c .+  egrid_c.*w .- xgridp_c
neg_ind_c = cons_c.<=zeros(nx,ne,nx)
ut_c = (cons_c.^(1-ssigma))./(1-ssigma)
function solve_age_c(currentState::modelState_y)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P_g       = currentState.P_b
    ut_g   = currentState.ut_b
    neg_ind_g   = currentState.neg_ind_b
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    #P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])

    # Calculate expected utility
    expected = zeros(nx,ne,nx)
    if (age < T)
        expected = repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx])
    end

    # Manipulate Matrices Budget constraint
    #xgrid_g  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
    #xgridp_g = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
    #egrid_g  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])

    # Compute consumption
    #cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
    # Compute Utility
    utility = ut_g .+ bbeta.*expected
    # Eliminate negatives
    #neg_ind = cons.<=zeros(nx,ne,nx)
    utility = utility.*(.!neg_ind_g).+neg_ind_g.*(-10.0^(5))
    VV = reshape(maximum(utility,dims = 3),(nx,ne)) #reshape is allowed

    return VV
end
function solve_VFI_c(T,ne,nx,P_c,ut_c,neg_ind_c,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
    for age = T:-1:1
        currentState = modelState_y(ne,nx,T,age,P_c,ut_c,neg_ind_c,ssigma,bbeta,V_tomorrow_cpu,w,r)
        tempV_cpu = solve_age_c(currentState)
        V_cpu[age,:,:] = tempV_cpu
        V_tomorrow_cpu = tempV_cpu
    end

    return V_cpu
end

function solve_age_y(currentState::modelState_yy)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P       = currentState.P_b
    ut_g   = currentState.ut_b
    neg_ind_g   = currentState.neg_ind_b
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
        expected = zeros(nx,1)
        if (age < T)
            #expected = repeat(sum(P[ie,:].*V,dims=2),outer = [1,1,nx])
            expected = reshape(sum(V.*reshape(P[ie,:],(1,ne)),dims=2),(nx,1))
        end

        # Manipulate Matrices Budget constraint
        #xgrid_g  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
        #xgridp_g = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
        #egrid_g  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])

        # Compute consumption
        #cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
        # Compute Utility
        utility = reshape(ut_g[:,ie,:],(nx,nx)) .+ bbeta.*expected
        # Eliminate negatives
        #neg_ind = cons.<=zeros(nx,ne,nx)
        utility = utility.*(.!reshape(neg_ind_g[:,ie,:],(nx,nx))).+reshape(neg_ind_g[:,ie,:],(nx,nx)).*(-10.0^(5))
        VV[:,ie] = reshape(maximum(utility,dims = 2),(nx,1)) #reshape is allowed
    end
    return VV
end

function solve_VFI_y(T,ne,nx,P,ut_c,neg_ind_c,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
    for age = T:-1:1
        currentState = modelState_yy(ne,nx,T,age,P,ut_c,neg_ind_c,ssigma,bbeta,V_tomorrow_cpu,w,r)
        tempV_cpu = solve_age_y(currentState)
        V_cpu[age,:,:] = tempV_cpu
        V_tomorrow_cpu = tempV_cpu
    end

    return V_cpu
end

function solve_age_d(ne,nx,T,age,P_g,ut_g,neg_ind_g,ssigma,bbeta,V,w,r)
    #=
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P_g       = currentState.P_b
    ut_g   = currentState.ut_b
    neg_ind_g   = currentState.neg_ind_b
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    =#
    # Manipulate Matrices
    #P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])

    # Calculate expected utility
    expected = zeros(nx,ne,nx)
    if (age < T)
        expected = repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx])
    end

    # Manipulate Matrices Budget constraint
    #xgrid_g  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
    #xgridp_g = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
    #egrid_g  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])

    # Compute consumption
    #cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
    # Compute Utility
    utility = ut_g .+ bbeta.*expected
    # Eliminate negatives
    #neg_ind = cons.<=zeros(nx,ne,nx)
    utility = utility.*(.!neg_ind_g).+neg_ind_g.*(-10.0^(5))
    VV = reshape(maximum(utility,dims = 3),(nx,ne)) #reshape is allowed

    return VV
end
function solve_VFI_d(T,ne,nx,P_c,ut_c,neg_ind_c,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
    for age = T:-1:1
        #currentState = modelState_y(ne,nx,T,age,P_c,ut_c,neg_ind_c,ssigma,bbeta,V_tomorrow_cpu,w,r)
        tempV_cpu = solve_age_d(ne,nx,T,age,P_c,ut_c,neg_ind_c,ssigma,bbeta,V_tomorrow_cpu,w,r)
        V_cpu[age,:,:] = tempV_cpu
        V_tomorrow_cpu = tempV_cpu
    end

    return V_cpu
end
# I culd take those out and create a structure that includes those guys
function solve_age(currentState::modelState_a)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P       = currentState.P
    xgrid   = currentState.xgrid
    egrid   = currentState.egrid
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])

    # Calculate expected utility
    expected = zeros(nx,ne,nx)
    if (age < T)
        expected = repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx])
    end

    # Manipulate Matrices Budget constraint
    xgrid_g  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
    xgridp_g = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
    egrid_g  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])

    # Compute consumption
    cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
    # Compute Utility
    utility = (cons.^(1-ssigma))/(1-ssigma) + bbeta.*expected
    # Eliminate negatives
    neg_ind = cons.<=zeros(nx,ne,nx)
    utility = utility.*(.!neg_ind)+neg_ind.*(-10.0^(5))
    VV = reshape(maximum(utility,dims = 3),(nx,ne)) #reshape is allowed

    return VV
end
function solve_VFI(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
    for age = T:-1:1
        currentState = modelState_a(ne,nx,T,age,P,xgrid,egrid,ssigma,bbeta,V_tomorrow_cpu,w,r)
        tempV_cpu = solve_age(currentState)
        V_cpu[age,:,:] = tempV_cpu
        V_tomorrow_cpu = tempV_cpu
    end

    return V_cpu
end


# only one staying should be V


# The original thread
function value_cpu(currentState::modelState_b)

  ind     = currentState.ind
  age     = currentState.age
  ne      = currentState.ne
  nx      = currentState.nx
  T       = currentState.T
  P       = currentState.P
  xgrid   = currentState.xgrid
  egrid   = currentState.egrid
  ssigma  = currentState.ssigma
  bbeta   = currentState.bbeta
  w       = currentState.w
  r       = currentState.r
  V       = currentState.V

  ix      = convert(Int, floor((ind-0.05)/ne))+1
  ie      = convert(Int, floor(mod(ind-0.05, ne))+1)

  VV      = -10.0^3
  ixpopt  = 0


    for ixp = 1:nx

      expected = 0.0
      if (age < T)
        for iep = 1:ne
          expected = expected + P[ie, iep]*V[ixp, iep]
        end
      end

      cons  = (1 + r)*xgrid[ix] + egrid[ie]*w - xgrid[ixp]

      utility = (cons^(1-ssigma))/(1-ssigma) + bbeta*expected # Change this for the point
      #utility = (cons^(1-ssigma))/(1-ssigma) .+ bbeta*expected
      if (cons <= 0)
        utility = -10.0^(5)
      end

      if (utility >= VV)
          # So I find the maximum over the ixp
        VV = utility # this is finding the maximum in a very funny way
        ixpopt = ixp
      end

      utility = 0.0
    end

    return VV

end

function VFI_cpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu,tempV_cpu)
    for age = T:-1:1

      Threads.@threads for ind = 1:(ne*nx)

        ix      = convert(Int, ceil(ind/ne))
        ie      = convert(Int, floor(mod(ind-0.05, ne))+1)

        currentState = modelState_b(ind,ne,nx,T,age,P,xgrid,egrid,ssigma,bbeta, V_tomorrow_cpu,w,r)
        tempV_cpu[ind] = value_cpu(currentState)
        V_cpu[age, ix, ie] = tempV_cpu[ind]
        V_tomorrow_cpu[ix, ie] = tempV_cpu[ind]
      end

      #finish = convert(Int, Dates.value(Dates.unix2datetime(time())- start))/1000;
      #print("Age: ", age, ". Time: ", finish, " seconds. \n")
    end
    return V_cpu
end
# GPU Alternative loading everything before hand
function solve_age_gpu(currentState::modelState_a)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P       = currentState.P
    xgrid   = currentState.xgrid
    egrid   = currentState.egrid
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])

    # Calculate expected utility
    expected = AFArray(zeros(nx,ne,nx))
    if (age < T)
        expected = AFArray(repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx]))
    end

    # Manipulate Matrices Budget constraint
    xgrid_g  = AFArray(permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2]))
    xgridp_g = AFArray(permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2]))
    egrid_g  = AFArray(repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx]))

    # Compute consumption
    cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
    # Compute Utility
    utility = (cons.^(1-ssigma))./(1-ssigma) + bbeta.*expected
    # Eliminate negatives
    neg_ind = cons.<=AFArray(zeros(nx,ne,nx))
    utility = utility.*(.!neg_ind)+neg_ind.*(-10.0^(5))
    VV = Array(reshape(maximum(utility,dims = 3),(nx,ne))) #reshape is allowed

    return VV
end

function solve_VFI_gpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_gpu,V_gpu)
    for age = T:-1:1
        currentState = modelState_a(ne,nx,T,age,P,xgrid,egrid,ssigma,bbeta,V_tomorrow_gpu,w,r)
        tempV_gpu = solve_age_gpu(currentState)
        V_gpu[age,:,:] = tempV_gpu
        V_tomorrow_gpu = tempV_gpu
    end

    return V_gpu
end
# Alternative take the matrices out
P_a = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
xgrid_a  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
xgridp_a = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
egrid_a  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])
cons_a = (1+r).*xgrid_a .+  egrid_a.*w .- xgridp_a
neg_ind_a = convert(Array{Bool},cons_a.<=zeros(nx,ne,nx)) # convert(Array{Int}, a)
ut_a = (cons_a.^(1-ssigma))./(1-ssigma)
function solve_age_gpu_a(currentState::modelState_x)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P_g       = currentState.P_a
    xgrid_g   = AFArray(currentState.xgrid_a)
    xgridp_g  = AFArray(currentState.xgridp_a)
    egrid_g   = AFArray(currentState.egrid_a)
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    #P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])

    # Calculate expected utility
    expected = AFArray(zeros(nx,ne,nx))
    if (age < T)
        expected = AFArray(repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx]))
    end

    # Manipulate Matrices Budget constraint
    #xgrid_g  = AFArray(permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2]))
    #xgridp_g = AFArray(permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2]))
    #egrid_g  = AFArray(repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx]))

    # Compute consumption
    cons = (1+r).*xgrid_g .+  egrid_g.*w .- xgridp_g
    # Compute Utility
    utility = (cons.^(1-ssigma))./(1-ssigma) + bbeta.*expected
    # Eliminate negatives
    neg_ind = cons.<=AFArray(zeros(nx,ne,nx))
    utility = utility.*(.!neg_ind).+neg_ind.*(-10.0^(5))
    VV = Array(reshape(maximum(utility,dims = 3),(nx,ne))) #reshape is allowed

    return VV
end

function solve_VFI_gpu_a(T,ne,nx,P_a,xgrid_a,xgridp_a,egrid_a,ssigma,bbeta,w,r,V_tomorrow_gpu_a,V_gpu_a)
    for age = T:-1:1
        currentState = modelState_x(ne,nx,T,age,P_a,xgrid_a,xgridp_a,egrid_a,ssigma,bbeta,V_tomorrow_gpu_a,w,r)
        tempV_gpu_a = solve_age_gpu_a(currentState)
        V_gpu_a[age,:,:] = tempV_gpu_a
        V_tomorrow_gpu_a = tempV_gpu_a
    end

    return V_gpu_a
end
function solve_age_gpu_b(currentState::modelState_zz)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P_g       = currentState.P_a
    ut_g      = AFArray(currentState.ut_a)
    neg_ind_g = AFArray(currentState.neg_ind_a)
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    #P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])

    # Calculate expected utility
    expected = AFArray(zeros(nx,ne,nx))
    if (age < T)
        expected = AFArray(repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx]))
    end

    # Manipulate Matrices Budget constraint
    #xgrid_g  = AFArray(permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2]))
    #xgridp_g = AFArray(permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2]))
    #egrid_g  = AFArray(repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx]))

    # Compute consumption
    #cons = (1+r).*xgrid_g .+  egrid_g.*w .- xgridp_g
    # Compute Utility
    utility = ut_g .+ bbeta.*expected
    # Eliminate negatives
    #neg_ind = cons.<=AFArray(zeros(nx,ne,nx))
    utility = utility.*(.!neg_ind_g).+neg_ind_g.*(-10.0^(5))
    VV = Array(reshape(maximum(utility,dims = 3),(nx,ne))) #reshape is allowed

    return VV
end

function solve_VFI_gpu_b(T,ne,nx,P_a,ut_a,neg_ind_a,ssigma,bbeta,w,r,V_tomorrow_gpu_b,V_gpu_b)
    for age = T:-1:1
        currentState = modelState_zz(ne,nx,T,age,P_a,ut_a,neg_ind_a,ssigma,bbeta,V_tomorrow_gpu_b,w,r)
        tempV_gpu_b = solve_age_gpu_b(currentState)
        V_gpu_b[age,:,:] = tempV_gpu_b
        V_tomorrow_gpu_b = tempV_gpu_b
    end

    return V_gpu_b
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
# CPU
# Original Thread Code
cpu_thread = @elapsed V_cpu_thread = VFI_cpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu,tempV_cpu)
# Vectorise ALL inside
cpu_vector = @elapsed V_cpu_solved = solve_VFI(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
# Vectorize all, Load utility
cpu_vector_c = @elapsed V_cpu_solved_c = solve_VFI_c(T,ne,nx,P_c,ut_c,neg_ind_c,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
# Vectorize only Savings
cpu_vector_y = @elapsed V_cpu_solved_y = solve_VFI_y(T,ne,nx,P,ut_c,neg_ind_c,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
# Vectorize all dont use structure
cpu_vector_d = @elapsed V_cpu_solved_d = solve_VFI_d(T,ne,nx,P_c,ut_c,neg_ind_c,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
## GPU
# Vectorize All, all inside (It will take triple time)
gpu_vector = @elapsed V_gpu_solved = solve_VFI_gpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_gpu,V_gpu)
# Vectorize All Load All (Takes too long)
gpu_vector_a = @elapsed V_gpu_solved_a = solve_VFI_gpu_a(T,ne,nx,P_a,xgrid_a,xgridp_a,egrid_a,ssigma,bbeta,w,r,V_tomorrow_gpu_a,V_gpu_a)
# Vectorize All Load Utility (Stack Overflow)
gpu_vector_b = @elapsed V_gpu_solved_b = solve_VFI_gpu_b(T,ne,nx,P_a,ut_a,neg_ind_a,ssigma,bbeta,w,r,V_tomorrow_gpu_b,V_gpu_b)
# Vectorize Savings
gpu_vector_c = @elapsed V_gpu_solved_c = solve_VFI_gpu_c(T,ne,nx,P,ut_a,neg_ind_a,ssigma,bbeta,w,r,V_tomorrow_gpu_b,V_gpu_b)

# Chech that results are the same
a_c1=V_cpu_thread.==V_cpu_solved
a_c2=V_cpu_thread.==V_cpu_solved_c
a_c3=V_cpu_thread.==V_gpu_solved
V_cpu_thread.==V_gpu_solved_a
V_cpu_thread.==V_gpu_solved_c
@profile solve_VFI_c(T,ne,nx,P_c,ut_c,neg_ind_c,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu)
