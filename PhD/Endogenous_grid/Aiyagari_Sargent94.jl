using Plots, Interpolations, LinearAlgebra, Parameters
#using Statistics, Compat, Expectations, NLsolve, Roots, Random, StatPlots # Distributions and Interpolations posible failure
plotly()
# functions

param = @with_kw (
        beta  = 0.96, # Discount Factor
        gamma = 2.0,  # Risk Aversion (inverse)
        alpha1 = 1/3,  # Capital Share
        delta = 0.08, # Depreciation Rate
        r = 0.036,
        w = 0.1,       # Wage

        # GRID Parameters

        na = 100, # more unknowns
        amax = 10,
        #a_temp = log.(range(exp(0.0), na, length=na)),
        a_temp = log.(range(exp(0.0), na, length=na)),
        agrid = (cumsum(a_temp)/sum(a_temp) * amax), # Asset Grid
        a_finer_temp  = log.(range(exp(0.0), 3*na, length=3*na)),
        agrid_finer   = (cumsum(a_finer_temp)/sum(a_finer_temp) * amax),
        b = agrid[1],

        # Iteration options
        dist_tol = 1e-5,
        maxit = 1000,
        #lb = beta+2, # Lower bound
        lb = 0.037, # Lower bound
        ub = 1.0/beta-1.0-0.0001)

function rouwenhorst(;
                    rho = 0.9,
                    p = (1.0+0.9)/2, #p = (1.0+rho)/2,
                    q = p,
                    sigma_e = 0.1, # Dont know
                    sigma_y = sqrt((sigma_e^2)/(1-rho^2)),
                    ny = 4, # NUmber of states
                    my = 0) # Who knows) # Number of states

        phy = sigma_y*sqrt(ny-1)
        ymax = phy
        ymin = -phy
        ygrid = range(ymin, ymax, length=ny)

        Piy2 = [[p 1.0-p];[1.0-q q]]
        Piyn1 = copy(Piy2)

        for jj = 1:(ny-2)
            num_rows = size(Piyn1,1)
            mat1     = zeros(num_rows+1, num_rows+1)
            mat2, mat3, mat4 = copy(mat1), copy(mat1), copy(mat1)

            mat1[1:end-1, 1:end-1]  = Piyn1
            mat2[1:end-1, 2:end]    = Piyn1
            mat3[2:end, 1:end-1]    = Piyn1
            mat4[2:end, 2:end]      = Piyn1

            Piyn1 = p*mat1 + (1-p)*mat2 + (1-q)*mat3 + q*mat4
            Piyn1[2:end-1, :] = Piyn1[2:end-1, :] / 2
        end
        Piy     = copy(Piyn1)
        Piy_aux = copy(Piy')
        vals = eigvals(Piy_aux)
        vecs = eigvecs(Piy_aux)
        todelete, ind_val  = findmin(abs.(vals.-1.0))
        Pyinv       = vecs[:, ind_val]/sum(vecs[:, ind_val])

        sum(Pyinv.>=0.0) == ny || throw(error("Negative elements in invariant distribution"))
        egrid   = exp.(ygrid + my*ones(ny))
        return (Pyinv=Pyinv, Piy=Piy ,egrid=egrid ,ny=ny, ygrid=ygrid)
        #return Pyinv
           # Persistence of the income process
end
#

function compute_aiyagari(mcm;)
    @unpack beta, delta, na, agrid, agrid_finer, gamma, alpha1 = mcm
    @unpack r, w, maxit, dist_tol, b = mcm
    #@unpack beta, delta, na, agrid, agrid_finer, gamma, alpha1 = param()
    #@unpack r, w, maxit, dist_tol, b = param()
    @unpack Pyinv, Piy, egrid, ny, ygrid = rouwenhorst()
    H = dot(egrid, Pyinv)
    #K = 0.0
    income_grid = egrid*w

    snodes      = [repeat(agrid_finer, ny) kron(ygrid, agrid_finer)]
    num_nodes   = size(snodes, 1)

    ### Utility Function and F.O.C. ###
    uprime(c) = c.^(-gamma)
    upinv(vv) = vv.^(-1.0/gamma)

    ### Production Function and F.O.C. ###
    prod_fun(K) = K^alpha1*H^(1.0-alpha1)                               # computes output
    mpk_fun(K)  = (alpha1*K^(alpha1-1.0)*H^(1.0-alpha1))                 # computes r
    mpl_fun(K)  = (1.0-alpha1)*K^alpha1*H^-alpha1                        # computes w
    inv_mpk(rate1) = H^(1.0-alpha1)*(alpha1/(rate1+delta))^(1.0/(1.0-alpha1))  # computes k

    # Initial Values
    r0_fixed = -0.1*delta + 0.9*(1.0/beta - 1.0)
    KK0   = inv_mpk(r0_fixed)
    ww0   = mpl_fun(KK0)
    income0 = egrid*ww0

    # Initialize consumption policy function
    cpol_mat = zeros(na, ny)
    for (i_z, z_v) in enumerate(income0)
        for (i_a, a_v) in enumerate(agrid)
            c_max = (1.0 + r0_fixed)*a_v + ww0*z_v-0.05
            cpol_mat[i_a, i_z] = c_max
        end
    end
    cpol_next = Array{Float64}(undef, na, ny)
    a_ast_mat = repeat(agrid, 1, ny)
    R = 1.0+r
    iter = 0
    dist = 10.0
    a_ast = Array{Float64}(undef, na, ny)
    while iter < maxit && dist > dist_tol
        iter +=1
        #old_pol = copy(cpol_mat)

        B = (beta * R * Piy * uprime(cpol_mat'))' # RHS of the Bellman Equation
        ctil  = upinv(B) # Solve for the ctil policy
        # Compute the endoegenous grid. Save it for piecewise linear approximation of the invariant distribution
        a_ast=(repeat(agrid, 1, ny) + ctil - repeat(income_grid', na, 1))/R # value of assets today that induces aprime tomorrow

        for (i_z, z_v) in enumerate(income_grid)
           for (i_a, a_v) in enumerate(agrid)
                if a_v<= a_ast[1, i_z] # case where you are constrained
                  cpol_next[i_a, i_z] = R*a_v + b + z_v
                elseif  a_v >= a_ast[end, i_z] # out of the range (to the right), linearly extrapolate
                  cpol_next[i_a,i_z] = ctil[end, i_z] + (a_v-a_ast[end, i_z])*(ctil[end, i_z] - ctil[end-1, i_z])/(a_ast[end, i_z]-a_ast[end-1, i_z])
                else # inside the range, linearly interpolate
                  ind1 = searchsortedfirst(a_ast[:,i_z], a_v)
                  ind  = ind1-1
                  cpol_next[i_a, i_z] = ctil[ind, i_z] + (a_v-a_ast[ind, i_z])*(ctil[ind1, i_z] - ctil[ind, i_z])/(a_ast[ind1, i_z]-a_ast[ind, i_z])
                end
           end
        end
        #dist = maximum(abs.(old_pol - cpol_mat))
        dist = norm(cpol_next - cpol_mat)
        cpol_mat .= cpol_next
    end
#=
    # Plots

    ap_egm = (1.0+ay.r) * repeat(ay.agrid, 1, ay.ny) + ay.w*repeat(ay.ϵgrid', ay.na, 1) - ay.cpol_mat
    fig, axis = subplots(1, 2)
    ax = axis[1]; ax[:plot](ay.agrid, ap_egm);      ax[:set_title]("Assets Policy at Equilibrium")
    ax = axis[2]; ax[:plot](ay.agrid, ay.cpol_mat); ax[:set_title]("Consumption Policy at Equilibrium")
    =#

    return (cpol_mat=cpol_mat, a_ast=a_ast)
    #return w
end
compute_aiyagari(param())
#=
function solveaiyagari(param;
    )
rouwenhorst()
H = dot(egrid, Pyinv)
K = 0

cpol_mat = zeros(param.na)
=#
