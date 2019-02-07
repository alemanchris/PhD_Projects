using InstantiateFromURL
activate_github("QuantEcon/QuantEconLecturePackages", tag = "v0.9.5");
using LinearAlgebra, Statistics, Compat, Interpolations
using Parameters, Plots, QuantEcon
using Polynomials
plotly()
Household = @with_kw (r = 0.03,
                      w = 1.0,
                      σ = 1.0,
                      gamma = 1.0,
                      β = 0.96,
                      delta = 0.05,
                      alpha1 = 0.33,
                      beta = 0.96,
                      #z_chain = MarkovChain([0.9 0.1; 0.1 0.9], [0.1; 1.0]),
                      z_chain = MarkovChain([0.7 0.3; 0.3 0.7], [0.1; 1.0]),
                      Pyinv = stationary_distributions(z_chain)[:, 1][1],
                      Piy = [0.7 0.3; 0.3 0.7],
                      #Piy = [0.9 0.1; 0.1 0.9],
                      z_vals = [0.1 1.0],
                      egrid =[z_vals[1];z_vals[2]],
                      a_min = 1e-10,
                      #a_min = -4,
                      b = a_min,
                      dist_tol = 1e-7,
                      #a_max = 20.0,
                      a_max = 20.0,
                      maxit = 2000,
                      a_size = 200,#200,
                      na = a_size,
                      a_size2 =3000,#250 #550
                      a_vals = range(a_min, a_max, length = a_size),
                      agrid_finer = range(a_min, a_max, length = a_size2),
                      agrid = a_vals,
                      z_size = length(z_chain.state_values),
                      ny = z_size,
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
function c_interp(Amat,cpo,val)
    c_aux = interp(Amat, cpo)
    return c_aux.(val)
end

function compute_invariant(mcm,a_ast)
    # For this particula problem a_ast will be the asset policy function
    # retrieving parameters
    @unpack na, agrid, agrid_finer= mcm
    @unpack maxit, dist_tol, b = mcm
    @unpack Pyinv, Piy, egrid, ny = mcm
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
    a_ast_itp = LinearInterpolation((agrid, egrid), a_ast) # next period's assets, current y
#BSpline(Linear())
    while iter<2500 && dist>1e-4
    #while dist>dist_tol
        iter += 1
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
        dist = norm(Λn_mat - Λnm1_mat)
        Λnm1_mat .= Λn_mat


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

        A_supply += sum_val_a + Λnm1_mat[1, i_y]*agrid_finer[1]
    end
    length_aux =length(agrid_finer)
    cdf_a = zeros(length_aux+1)

    for i_a=1:length_aux
        cdf_a[i_a] = sum(Λnm1_mat[i_a, :])
    end
    cdf_a[length_aux+1] = cdf_a[length_aux] #repeat the last value with the same so that the pdf of the last is zero
    # GET THE PDF
    pdf_a = diff(cdf_a)
    #mean_K2_aux=0
    #for i_y = 1:ny
    #    mean_K2_aux += dot(agrid_finer,Λnm1_mat[:,i_y])
    #end
    # Amother way to get the PDF
    index_pdf = pdf_a.==0
    dert2_aux = sum(index_pdf)
    pdf_d = pdf_a[1:end-dert2_aux]
    pdf_d = pdf_d./(sum(pdf_d))
    agrdi_finer3 = range(0,20,length=length(pdf_d))
    A_supply3 = dot(agrdi_finer3,pdf_d)
    return (A_supply3 = A_supply3, A_supply = A_supply, pdf_a=pdf_a, iter = iter, cdf_a=cdf_a, pdf_d=pdf_d, agrid_finer=agrid_finer)
end

function m_invariant(mcm,a_STAR)
    @unpack na, agrid, agrid_finer= mcm
    @unpack maxit, dist_tol, b = mcm
    @unpack Pyinv, Piy, egrid, ny,a_min,a_max = mcm
    Lambda0 = zeros(na*2,ny)
    for j = 1:ny
        for k = 1:na*2
            Lambda0(k,j) = (agrid_finer[k]-a_min)/(a_max-a_min)*Pyinv[j]
        end
    end
    # Make it continous
    #l0 = @(a,lambda) interp1(agrid2,lambda,a,'linear'); %  I can condition de densities
    #l0(agrid_finer,lambda,a)
    # update distribution for every pair
    #crit1 = 10^(-8);

    lambda = Lambda0
    l1 = zeros(na*2,ny);

    for i=1:ny
        for j=1:na*2
            #aux = interp(agrid,a_STAR[:,i],agrid_finer[j])
            #aux2 = interp1(kap,condecis(:,i),agrid_finer[j])
            aprime2[j,i] = maximum([minimum([i_STAR(agrid,a_STAR[:,i],agrid_finer[j]),a_max]),a_min])
            #condecis2(j,i) = max(aux2,0);
        end
    end
    #maxiter = 10
    #Inverse mapping
    #invaprime(a,i) =interp(a_STAR[:,i],agrid_finer,a)
    iter3 = 0
    dist1 = 1
    #range  = [a_min,amax]
    #while iter<2000 && dist1>1e-4
    while iter<1000 && dist1>1e-4

        iter3 = iter3+1
            for k = 1:na*2
                for j = 1:ny
                    for i = 1:ny
                        #%dert= invaprime(a,s,r,agrid2',condecis2,agrid2(k),i,wage),X0)
                        #asol = fzero(invaprime(a,s,r,agrid2',condecis2,agrid2(k),i,wage))
                        #asol = invaprime(agrid_finer[k],i)
                        asol = invaprime(agrid_finer[k],a_STAR[:,i])
                        asol2= maximum([minimum([asol,a_max]),a_min])
                        l1[k,j] = l1[k,j]+prob[i,j]*l0(agrid_finer,lambda[:,i],asol2)
                    end
                end
            end


        dist1 = maximum(maximum(abs.(lambda-l1)))
        #dist1 = norm(abs.(lambda-l1)
        lambda .= l1
    end
    A_supply2 = 0
    # This is calculating an integral manually
    for i_y = 1:ny
        sum_val_a = 0.0
        for i_a = 1:(length(agrid_finer)-1)
            anp1 = agrid_finer[i_a+1]
            an   = agrid_finer[i_a]
            sum_val_a += 0.5*(lambda[i_a+1, i_y] - lambda[i_a, i_y]) * (anp1 + an)
        end

        A_supply2 += sum_val_a + lambda[1, i_y]*agrid_finer[1]
    end

    length_aux =length(agrid_finer)
    cdf_b = zeros(length_aux+1)

    for i_a=1:length_aux
        cdf_b[i_a] = sum(lambda[i_a, :])
    end
    cdf_b[length_aux+1] = cdf_b[length_aux] #repeat the last value with the same so that the pdf of the last is zero
    # GET THE PDF
    pdf_b = diff(cdf_b)

    return (A_supply3 = A_supply3, pdf_d=pdf_d, iter3= iter3 )
end

function lo(agrid_finer,lambda,a)
    dert_aux = interp(agrid_finer,lambda)
    return dert_aux(a)
end

function i_STAR(agrid,a_STAR,agrid_finer)
    dert_aux2 = interp(agrid,a_STAR)
    return dert_aux2(agrid_finer)
end
function invaprime(a,a_STAR,agrid_finer)
    dert_aux3 =interp(a_STAR,agrid_finer)
    return dert_aux3(a)
end
