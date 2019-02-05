using InstantiateFromURL
activate_github("QuantEcon/QuantEconLecturePackages", tag = "v0.9.5");
using LinearAlgebra, Statistics, Compat
using Parameters, Plots, QuantEcon
#gr(fmt = :png);
n = 7
rho = 0.2
sigmaint = 0.4           # intermediate value to calculate sigma
sigma    = sigmaint*sqrt(1-rho^2) # standard deviation of error_t
mc = tauchen(n, rho, sigma)
Pyinv = stationary_distributions(mc)[1]
invariant = stationary_distributions(mc)
dert = log.(mc[1,:])
mc[1,:]
trans = mc.p

##### Arguments
#=
- `N::Integer`: Number of points in markov process
- `ρ::Real` : Persistence parameter in AR(1) process
- `σ::Real` : Standard deviation of random component of AR(1) process
- `μ::Real(0.0)` : Mean of AR(1) process
- `n_std::Integer(3)` : The number of standard deviations to each side the process
  should span
##### Returns
- `mc::MarkovChain{Float64}` : Markov chain holding the state values and transition matrix
"""
function tauchen(N::Integer, ρ::Real, σ::Real, μ::Real=0.0, n_std::Integer=3)
=#
