using InstantiateFromURL
activate_github("QuantEcon/QuantEconLecturePackages", tag = "v0.9.5");
using LinearAlgebra, Statistics, Compat
using Parameters, Plots, QuantEcon
Household = @with_kw (r = 0.01,
                      w = 1.0,
                      σ = 1.0,
                      β = 0.96,
                      z_chain = MarkovChain([0.9 0.1; 0.1 0.9], [0.1; 1.0]),
                      Piy = [0.9 0.1; 0.1 0.9],
                      a_min = 1e-10,
                      a_max = 18.0,
                      a_size = 200,
                      a_vals = range(a_min, a_max, length = a_size),
                      z_size = length(z_chain.state_values),
                      n = a_size * z_size,
                      s_vals = gridmake(a_vals, z_chain.state_values),
                      s_i_vals = gridmake(1:a_size, 1:z_size),
                      u = σ == 1 ? x -> log(x) : x -> (x^(1 - σ) - 1) / (1 - σ),
                      R = setup_R!(fill(-Inf, n, a_size), a_vals, s_vals, r, w, u),
                      # -Inf is the utility of dying (0 consumption)
                      Q = setup_Q!(zeros(n, a_size, n), s_i_vals, z_chain))

function setup_Q!(Q, s_i_vals, z_chain)
    for next_s_i in 1:size(Q, 3)
        for a_i in 1:size(Q, 2)
            for s_i in 1:size(Q, 1)
                z_i = s_i_vals[s_i, 2]
                next_z_i = s_i_vals[next_s_i, 2]
                next_a_i = s_i_vals[next_s_i, 1]
                if next_a_i == a_i
                    Q[s_i, a_i, next_s_i] = z_chain.p[z_i, next_z_i]
                end
            end
        end
    end
    return Q
end

function setup_R!(R, a_vals, s_vals, r, w, u)
    for new_a_i in 1:size(R, 2)
        a_new = a_vals[new_a_i]
        for s_i in 1:size(R, 1)
            a = s_vals[s_i, 1]
            z = s_vals[s_i, 2]
            c = w * z + (1 + r) * a - a_new
            if c > 0
                R[s_i, new_a_i] = u(c)
            end
        end
    end
    return R
end
const A = 1
const N = 1
const α = 0.33
const β = 0.96
const δ = 0.05

function r_to_w(r)
    return A * (1 - α) * (A * α / (r + δ)) ^ (α / (1 - α))
end

function prices_to_capital_stock(am, r)

    # Set up problem
    w = r_to_w(r)
    @unpack a_vals, s_vals, u = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)

    aiyagari_ddp = DiscreteDP(am.R, am.Q, am.β)

    # Compute the optimal policy
    results = solve(aiyagari_ddp, PFI)

    # Compute the stationary distribution
    stationary_probs = stationary_distributions(results.mc)[:, 1][1]

    # Return K
    return dot(am.s_vals[:, 1], stationary_probs)
end
