r0_fixed = -0.1*ay.δ + 0.9*(1.0/ay.β - 1.0)

""" Initialize consumption policy """
KK0   = inv_mpk(ay, r0_fixed)
ww0   = mpl_fun(ay, KK0)
zgrid = ay.ϵgrid*ww0

cpol_init = zeros(ay.na, ay.ny)

for (i_z, z_v) in enumerate(zgrid)
    for (i_a, a_v) in enumerate(ay.agrid)
        c_max = (1.0 + r0_fixed)*a_v + ww0*z_v-0.05
        cpol_init[i_a, i_z] = c_max
    end
end

copyto!(ay.cpol_mat, cpol_init)

#copyto!(ay.r, 0.031)
ay.r = 0.031
ay.K = inv_mpk(ay, ay.r)
ay.w = mpl_fun(ay, ay.K)

copyto!(ay.cpol_mat, cpol_init)

""" Solve using EndoGrid """
iterate_pol!(ay)


""" Compute invariant distribution given r0 """
Λ_invariant = zeros(length(ay.agrid_finer), ay.ny)

for i_y=1:ay.ny
    for (i_a, a_v) in enumerate(ay.agrid_finer)
        Λ_invariant[i_a, i_y] = (a_v - ay.agrid[1]) / (ay.agrid[end] - ay.agrid[1]) * ay.Πyinv[i_y]
    end
end

compute_invariant_pwlinear!(ay, Λ_invariant)
