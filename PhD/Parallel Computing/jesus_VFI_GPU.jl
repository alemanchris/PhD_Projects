using Distributions
using ArrayFire
#--------------------------------#
#         Initialization         #
#--------------------------------#


# Grid for x
nx  = 12 #1500
xmin = 0.1
xmax = 4.0

# Grid for e: parameters for Tauchen
ne  = 4 #15
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
V_gpu          = AFArray(zeros(T, nx, ne))
V_tomorrow_cpu = zeros(nx, ne)
V_tomorrow_gpu = AFArray(zeros(nx, ne))

# Initialize value function as a shared array
tempV_cpu = zeros(ne*nx)
tempV_gpu = AFArray(zeros(ne*nx))
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
  V::AFArray{Float64,2}
  w::Float64
  r::Float64
end

function value_gpu(currentState::modelState_b)

 ind     = currentState.ind
 age     = currentState.age
 ne      = currentState.ne
 nx      = currentState.nx
 T       = currentState.T
 P       = AFArray(currentState.P)
 xgrid   = AFArray(currentState.xgrid)
 egrid   = AFArray(currentState.egrid)
 ssigma  = currentState.ssigma
 bbeta   = currentState.bbeta
 w       = currentState.w
 r       = currentState.r
 V       = currentState.V

 #ix      = convert(Int, floor((ind-0.05)/ne))+1
 #ie      = convert(Int, floor(mod(ind-0.05, ne))+1)

 ix      = floor(Int,(ind-0.05)/ne)+1
 ie      = floor(Int,mod(ind-0.05, ne)+1)
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
     utility = (cons^(1-ssigma))/(1-ssigma) + bbeta*expected
     if (cons <= 0)
       utility = -10.0^(5)
     end

     if (utility >= VV)
       VV = utility
       ixpopt = ixp
     end
     utility = 0.0
   end
   return VV
end
#
function de_debug(T)
    dert  = 1
    for ind = 1:(ne*nx)

        ix      = convert(Int, ceil(ind/ne))
        ie      = convert(Int, floor(mod(ind-0.05, ne))+1)

        currentState = modelState_b(ind,ne,nx,T,T,P,xgrid,egrid,ssigma,bbeta, V_tomorrow_gpu,w,r)
        tempV_cpu[ind] = value_gpu(currentState)
        V_gpu[T, ix, ie] = tempV_gpu[ind]
        V_tomorrow_gpu[ix, ie] = tempV_gpu[ind]

    end
    dert2 = 3
    dert3 = 4
    return V_gpu
end
V_small = de_debug(T)
dert = 1
#@Juno.enter de_debug(T)
#@Juno.enter de_debug(T)
dert_15 = [1 1 1 1 1]
dert_16 = Vector{Float64,[1.;1.;1.;1.;1.]}
