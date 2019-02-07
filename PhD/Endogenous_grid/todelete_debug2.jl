# Create an Instance of a Household
am = Household(a_max = 20.0, r = 0.03, w = 0.956)
# Do manually
r = 0.03
w = r_to_w(r)
@unpack a_vals, s_vals, u, Piy,z_vals,z_chain = am
setup_R!(am.R, a_vals, s_vals, r, w, u)

aiyagari_ddp = DiscreteDP(am.R, am.Q, am.β)

# Compute the optimal policy
results = solve(aiyagari_ddp, PFI)
# Obtain the values of the policy functions
@unpack z_size, a_size, n = am
a_star = reshape([a_vals[results.sigma[s_i]] for s_i in 1:n], a_size, z_size)
a_star_idx = reshape([results.sigma[s_i] for s_i in 1:n], a_size, z_size)

# Debug
@unpack na, agrid, agrid_finer= Household(a_max = 20.0, r = 0.03)
@unpack maxit, dist_tol, b = Household(a_max = 20.0, r = 0.03)
@unpack Pyinv, Piy, egrid, ny= Household(a_max = 20.0, r = 0.03)
#@unpack cpol_mat, a_ast, apol_egm = compute_aiyagari(mcm)
# Initial Values
Λ_invariant_init = zeros(length(agrid_finer), ny)
#val_aux = 1/length(agrid_finer)*ny
#Λ_invariant_init = zeros(length(agrid_finer), ny).+val_aux
#
for i_y=1:ny
    for (i_a, a_v) in enumerate(agrid_finer)
        Λ_invariant_init[i_a, i_y] = (a_v - agrid[1]) / (agrid[end] - agrid[1]) * Pyinv[i_y]
    end
end
#
# Initialize Matrices
length_aux = length(agrid_finer)
Λn_mat   = zeros(length(agrid_finer), ny) # Distribution
Λnm1_mat = copy(Λ_invariant_init) # New distribution
# Start Iteration
iter = 0
dist = 10.0
a_ast_itp = LinearInterpolation((agrid, egrid), a_star) # next period's assets, current y

while iter<200 && dist>1e-4
#while dist>dist_tol
    global iter += 1
    for i_y = 1:ny # next period
        for (i_a, a_v) in enumerate(agrid_finer) # next period
            vec_temp = zeros(ny)
            for (i_y0, y0_v) in enumerate(egrid) # last period y
                aux_ast = a_ast_itp(a_v, y0_v)
                aval  = minimum([maximum([aux_ast, agrid[1]]), agrid[end]]) # today's assets (endogenous grid)
                aux_ast_pos = searchsortedfirst(agrid_finer, aval)
                ind_r = minimum([maximum([aux_ast_pos, 2]), length_aux])
                #
                Λval = Λnm1_mat[ind_r-1, i_y0] + (Λnm1_mat[ind_r, i_y0]- Λnm1_mat[ind_r-1, i_y0]) / (agrid_finer[ind_r] - agrid_finer[ind_r-1]) * (aval - agrid_finer[ind_r-1])
                vec_temp[i_y0] = Λval
            end

            Λn_mat[i_a, i_y] = dot(Piy[:, i_y], vec_temp) #the distribution
        end
    end

    #dist = maximum(abs.(Λn_mat - Λnm1_mat))
    #copyto!(Λnm1_mat, Λn_mat)
    global dist = norm(Λn_mat - Λnm1_mat)
    global Λnm1_mat .= Λn_mat


    # if iter%200 == 0.0
    #     @printf("Iteration # %d on distribution with distance %.3g\n", iter, dist)
    # end
end
A_supply = 0
# This is calculating an integral manually
for i_y = 1:ny
    sum_val_a = 0.0
    for i_a = 1:(length(agrid_finer)-1)
        anp1 = agrid_finer[i_a+1]
        an   = agrid_finer[i_a]
        sum_val_a += 0.5*(Λnm1_mat[i_a+1, i_y] - Λnm1_mat[i_a, i_y]) * (anp1 + an)
    end

    global A_supply += sum_val_a + Λnm1_mat[1, i_y]*agrid_finer[1]
end
