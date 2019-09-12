#--------------------------------#
#         House-keeping          #
#--------------------------------#

using Distributed
using Distributions
using ArrayFire

#--------------------------------#
#         Initialization         #
#--------------------------------#


# Grid for x
nx  = 1500
xmin = 0.1
xmax = 4.0

# Grid for e: parameters for Tauchen
ne  = 15
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
V_tomorrow_cpu = zeros(nx, ne)
V_tomorrow_gpu = zeros(nx, ne)

# Initialize value function as a shared array
tempV_cpu = zeros(ne*nx)
tempV_gpu = zeros(ne*nx)
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



#--------------------------------#
#     Structure and function     #
#--------------------------------#

# Data structure of state and exogenous variables
struct modelState_a
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
# Function that computes value function, given vector of state variables
function value_cpu(currentState::modelState_a)

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
function value_gpu(currentState::modelState_b)

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
 #ie = convert(Int, floor((ind-0.05)/nx))+1
 #ix      = convert(Int, floor(mod(ind-0.05, nx))+1)
 #ix      = floor(Int,(ind-0.05)/ne)+1
 #ie      = floor(Int,mod(ind-0.05, ne)+1)
 #VV      = -10.0^3
 ixpopt  = 0
# If I would foguet about looping over the ix

   #for ixp = 1:nx
 P_a = transpose(P)
     #expected = AFArray(zeros(1,nx))
     expected = zeros(1,nx)
     if (age < T)
         #
       for iep = 1:ne
         #expected = expected + P[ie, iep].*V[ixp, iep]
         expected = expected .+ transpose(P[ie, iep].*V[:, iep])
       end
       #expected = transpose(expected)
       #
       #expected = AFArray(reshape(sum(V.*reshape(P[ie,:],(1,ne)),dims=2),(1,nx))) # (nx,2) (nx,1) (1,nx)
       #expected = sum(V.*P[ie,:],dims=2)' # (nx,2) (nx,1) (1,nx)
     end

     #cons  = (1 + r)*xgrid[ix] + egrid[ie]*w - xgrid[ixp]
     #cons  = (1 + r).*xgrid[ix] .+ egrid[ie].*w .- AFArray(xgrid)
     #cons  = AFArray((1 + r).*xgrid[ix] .+ egrid[ie].*w .- xgrid)
     cons  = (1 + r).*xgrid[ix] .+ egrid[ie].*w .- xgrid
     utility = (cons.^(1-ssigma))./(1-ssigma) .+ bbeta.*expected
     neg_ind = cons.<=zeros(1,nx)
     utility = utility.*(.!neg_ind).+(neg_ind.*(-10.0^(5)))
     VV_aux = maximum(utility,dims=2)
     VV = VV_aux[1]
     #=
     if (cons <= 0)
       utility = -10.0^(5)
     end

     if (utility >= VV)
       VV = utility
       ixpopt = ixp
     end

     utility = 0.0
   end
   =#
   return VV

end
#=
function de_debug(T)

    for ind = 1:(ne*nx)

        ix      = convert(Int, ceil(ind/ne))
        ie      = convert(Int, floor(mod(ind-0.05, ne))+1)

        currentState = modelState_b(ind,ne,nx,T,T,P,xgrid,egrid,ssigma,bbeta, V_tomorrow_gpu,w,r)
        tempV_cpu[ind] = value_gpu(currentState)
        V_gpu[T, ix, ie] = tempV_gpu[ind]
        V_tomorrow_gpu[ix, ie] = tempV_gpu[ind]
        return V_gpu
    end
end
@Juno.enter de_debug(T)
=#
#--------------------------------#
#     Life-cycle computation     #
#--------------------------------#

print(" \n")
print("Life cycle computation: \n")
print(" \n")

#start = Dates.unix2datetime(time())
function VFI_cpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu,tempV_cpu)
    for age = T:-1:1

      Threads.@threads for ind = 1:(ne*nx)

        ix      = convert(Int, ceil(ind/ne))
        ie      = convert(Int, floor(mod(ind-0.05, ne))+1)

        currentState = modelState_a(ind,ne,nx,T,age,P,xgrid,egrid,ssigma,bbeta, V_tomorrow_cpu,w,r)
        tempV_cpu[ind] = value_cpu(currentState)
        V_cpu[age, ix, ie] = tempV_cpu[ind]
        V_tomorrow_cpu[ix, ie] = tempV_cpu[ind]
      end

      #finish = convert(Int, Dates.value(Dates.unix2datetime(time())- start))/1000;
      #print("Age: ", age, ". Time: ", finish, " seconds. \n")
    end
    return V_cpu
end

function VFI_gpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_gpu,V_gpu,tempV_gpu)
    for age = T:-1:1

      for ind = 1:(ne*nx)

        ix      = convert(Int, ceil(ind/ne))
        ie      = convert(Int, floor(mod(ind-0.05, ne))+1)
        #ie = convert(Int, floor((ind-0.05)/nx))+1
        #ix      = convert(Int, floor(mod(ind-0.05, nx))+1)
        currentState = modelState_b(ind,ne,nx,T,age,P,xgrid,egrid,ssigma,bbeta, V_tomorrow_gpu,w,r)
        tempV_gpu[ind] = value_gpu(currentState)
        V_gpu[age, ix, ie] = tempV_gpu[ind]
        V_tomorrow_gpu[ix, ie] = tempV_gpu[ind]
      end
      #=
      for ind = 1:(ne*nx)

        ix      = convert(Int, ceil(ind/ne));
        ie      = convert(Int, floor(mod(ind-0.05, ne))+1)

        V_gpu[age, ix, ie] = tempV_gpu[ind]
        V_tomorrow_gpu[ix, ie] = tempV_gpu[ind]
      end
      =#
      #finish = convert(Int, Dates.value(Dates.unix2datetime(time())- start))/1000;
      #print("Age: ", age, ". Time: ", finish, " seconds. \n")
    end
    return V_gpu
end


time_CPU = @elapsed V_cpu2 = VFI_cpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_cpu,V_cpu,tempV_cpu)
time_GPU = @elapsed V_gpu  = VFI_gpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_gpu,V_gpu,tempV_gpu)

V_cpu2 .== V_gpu
#print("\n")
#finish = convert(Int, Dates.value(Dates.unix2datetime(time())- start))/1000;
#print("TOTAL ELAPSED TIME: ", finish, " seconds. \n")

#---------------------#
#     Some checks     #
#---------------------#
#=
print(" \n")
print(" - - - - - - - - - - - - - - - - - - - - - \n")
print(" \n")
print("The first entries of the value function: \n")
print(" \n")

# I print the first entries of the value function, to check
for i = 1:3
  print(round(V[1, 1, i], digits=5), "\n")
end
=#
