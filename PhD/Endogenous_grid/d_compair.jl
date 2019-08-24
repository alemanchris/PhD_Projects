#=
 The following code uses a variaty of methods to approximate the
 invariant distribution of kapital supply.
 The methods used are the following:
    1. Discretization taken from Sargent's QuantEcon.
    2. Manual Discretization
    3. Simulation
    4. Piecewise Linear interpolation A
    5. Piecewise Linear interpolation B
    6. Collocation
=#
# 0. Generate an Instance of household
am = Household(β = β, a_max = 20.0)
r = 0.03
# 1. Discretization: Bechmark example taken from Sarget's QuantEcon
    # Kap_A = Mean Capital
    # SSdist_A = Invariant distribution
@unpack Kap_A, SSdist_A = dist_sargent_A(am, r)

# 2. Discretization: Manually
@unpack Kap_B, SSdist_B = dist_sargent_B(am, r)

# 3. Simulation
@unpack Kap_C, SSdist_C = simulation_C(am, r)

# 4. Piecewise Linear-Interpolation A. Notes by J.Violante
