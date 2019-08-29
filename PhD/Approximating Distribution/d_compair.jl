#=
 The following code uses a variaty of methods to Compute the
 invariant distribution of capital supply.
 The methods used are the following:
    1. Eigen Vector Method taken from Sargent's QuantEcon (Benchmark)
    2. Manual Eigen Vector Method (Trying to replicate sargent)
    3. Monte Carlo Simulation
    4. Discretization of the CDF A (approximate the CDF)
    5. Discretization of the PDF B (approximate the PDF)
    6. Collocation
=#
# 0. Generate an Instance of household
am = Household(β = β, a_max = 20.0)
r = 0.03
# 1. Eigen Vector Method:  Bechmark example taken from Sarget's QuantEcon
    # Kap_A = Mean Capital
    # SSdist_A = Invariant distribution
@unpack Kap_A, SSdist_A = dist_sargent_A(am, r)

# 2. Eigen Vector Method: Manually
@unpack Kap_B, SSdist_B = dist_sargent_B(am, r)

# 3. Simulation
@unpack Kap_C, SSdist_C = simulation_C(am, r)

# 4. Discretization of the CDF: A. Notes by J.Violante
@unpack Kap_E, SSdist_E = CDF_disc_E(am, r)

# 5. Discretization of the PDF: B. Notes by J.Violante
@unpack Kap_F, SSdist_F = PDF_disc_F(am, r)

# 6. Piecewise Linear-Interpolation (Rios Rull 1997) (Still to be corrected)
@unpack Kap_F, SSdist_F = PWLinear_D(am, r)

# 7. Collocation Chebyschev as in Windberry (Still to be done)
