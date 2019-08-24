using Plots, Interpolations, LinearAlgebra, Parameters
#using Statistics, Compat, Expectations, NLsolve, Roots, Random, StatPlots # Distributions and Interpolations posible failure
using InstantiateFromURL
activate_github("QuantEcon/QuantEconLecturePackages", tag = "v0.9.5");
using Statistics, QuantEcon,Compat , Random

plotly()
# functions

param = @with_kw (
        #=
        beta  = 0.97, # Discount Factor
        gamma = 3,  # Risk Aversion (inverse)
        alpha1 = 0.33,  # Capital Share
        delta = 0.05, # Depreciation Rate
        r = 0.036,
        w = (1-alpha1)*((alpha1/(r+delta))^alpha1)^(1/(1-alpha1)),       # Wage
        #w = 0.2,
        egrid = 0.156, #natural borroeing limit usually before ad/hoc constraint
        b_cons = 0.5,
        phi = minimum([b_cons,w*egrid[1]/r]),
        # GRID Parameters

         # more unknowns
        amin = -phi,
        amax = 20,
        inckap = 0.1,
        agrid = range(amin,step=inckap,stop=amax),
        na = length(agrid),
        agrid_finer = range(amin,step=0.05,stop=amax),
        b = -phi,
        =#
        #=
        beta  = 0.96, # Discount Factor
        gamma = 2.0,  # Risk Aversion (inverse)
        alpha1 = 1/3,  # Capital Share
        delta = 0.08, # Depreciation Rate
        r = 0.036,
        w = 0.1,       # Wage

        # GRID Parameters

         # more unknowns
        #
        na = 100,
        amax = 10,#10
        a_temp = log.(range(exp(0.0), na, length=na)),
        agrid = range(0, amax, length=na),
        #agrid = (cumsum(a_temp)/sum(a_temp) * amax), # Asset Grid
        a_finer_temp  = log.(range(exp(0.0), 2na, length=2*na)),
        agrid_finer  = range(0, amax, length=3*na),
        #agrid_finer   = (cumsum(a_finer_temp)/sum(a_finer_temp) * amax),
        #agrid_finer = agrid,

        b = agrid[1],
        =#
        beta  = 0.96, # Discount Factor
        gamma = 3.0,  # Risk Aversion (inverse)
        alpha1 = 0.36,# Capital Share
        delta = 0.08, # Depreciation Rate
        r = 0.036,
        #w = (1-alpha1)*((alpha1/(r+delta))^alpha1)^(1/(1-alpha1)),
        w = 0.2,       # Wage
        phy = minimum([0,w*0.301194211912202/r]),

        # GRID Parameters

         # more unknowns
        #
        na = 81,
        amax = 16,#16
        #a_temp = log.(range(exp(0.0), na, length=na)),
        agrid = range(-phy, amax, length=na),
        #agrid = (cumsum(a_temp)/sum(a_temp) * amax), # Asset Grid
        #a_finer_temp  = log.(range(exp(0.0), 2na, length=2*na)),
        agrid_finer  = range(-phy, amax, length=2*na),
        #agrid_finer   = (cumsum(a_finer_temp)/sum(a_finer_temp) * amax),
        #agrid_finer = agrid,

        b = agrid[1],

        # Iteration options
        dist_tol = 1e-7,
        maxit = 1000,
        #lb = beta+2, # Lower bound
        lb = 0.037, # Lower bound
        ub = 1.0/beta-1.0-0.0001)

#
function markov_aprox(;
                     ny = 7,
                     rho = 0.2,
                     sigmaint = 0.4,
                     sigma    = sigmaint*sqrt(1-rho^2),
                     mc = tauchen(ny, rho, sigma), # Generates same process as markovappr in matlab Sargent
                     Pyinv = stationary_distributions(mc)[1], #[1] is the method
                     Piy = mc.p,
                     ygrid = mc.state_values,
                     egrid = exp.(ygrid))

    return (Pyinv=Pyinv, Piy=Piy ,egrid=egrid ,ny=ny, ygrid=ygrid, mc=mc)
end



function compute_aiyagari(mcm;
                         plots=0)
    @unpack beta, delta, na, agrid, agrid_finer, gamma, alpha1 = mcm
    @unpack r, w, maxit, dist_tol, b = mcm
    #@unpack Pyinv, Piy, egrid, ny, ygrid = rouwenhorst()
    @unpack Pyinv, Piy, egrid, ny, ygrid = markov_aprox()
    # Only prob and s=egrid and Pyinv egrid is exp of ygrid=logs

    ########
    H = dot(egrid, Pyinv)
    K = 0.0
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
    #r0_fixed = 0.02
    KK0   = inv_mpk(r0_fixed)
    ww0   = mpl_fun(KK0)
    income0 = egrid*ww0

    # Initialize consumption policy function
    cpol_mat = zeros(na, ny)
    #
    for (i_z, z_v) in enumerate(income0)
        for (i_a, a_v) in enumerate(agrid)
            c_max = (1.0 + r0_fixed)*a_v + ww0*z_v-0.05
            #c_max = (1.0 + r0_fixed)*a_v + ww0*z_v+20000000000
            cpol_mat[i_a, i_z] = c_max
        end
    end
    #
    cpol_next = Array{Float64}(undef, na, ny)
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
                  # Why manually interpolate?
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
    # Policy function assets
    apol_egm = (1.0+r) * repeat(agrid, 1, ny) + w*repeat(egrid', na, 1)*H - cpol_mat

    # Plots
    if plots==1
        plot(agrid, apol_egm)
        plot(agrid, cpol_mat)
    #fig, axis = subplots(1, 2)
    #ax = axis[1]; ax[:plot](ay.agrid, ap_egm);      ax[:set_title]("Assets Policy at Equilibrium")
    #ax = axis[2]; ax[:plot](ay.agrid, ay.cpol_mat); ax[:set_title]("Consumption Policy at Equilibrium")
    end
    return (cpol_mat=cpol_mat, a_ast=a_ast, apol_egm=apol_egm)
    # cpol_mat, consumption policy function
    # a_ast, asset grid, found endogenously
    # apol_egm, asset policy function
    #return w
end

function compute_invariant(mcm;)
    # retrieving parameters
    @unpack beta, delta, na, agrid, agrid_finer, gamma, alpha1 = mcm
    @unpack r, w, maxit, dist_tol, b = mcm
    #@unpack Pyinv, Piy, egrid, ny, ygrid = rouwenhorst()
    @unpack Pyinv, Piy, egrid, ny, ygrid, mc= markov_aprox()

    income_grid = egrid*w
    @unpack cpol_mat, a_ast, apol_egm = compute_aiyagari(mcm,plots=1)
    # Initial Values
    H = dot(egrid, Pyinv)
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
    a_ast_itp = LinearInterpolation((agrid, ygrid), a_ast) # next period's assets, current y

    while iter<200 && dist>1e-4
    #while dist>dist_tol
        iter += 1
        for i_y = 1:ny # next period
            for (i_a, a_v) in enumerate(agrid_finer) # next period
                vec_temp = zeros(ny)
                for (i_y0, y0_v) in enumerate(ygrid) # last period y
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
        # Discretization Sanrgent
        # Getting the indices
        ind_asg = Array{Int64}(undef,na,ny)
        for i_z = 1:ny
            for i_a = 1:na
                  useles = findmin(abs.(agrid.-apol_egm[i_a,i_z]))
                  #ind_aux = minimum([searchsortedfirst(agrid, apol_egm[i_a,i_z]),na])
                  ind_asg[i_a,i_z] = useles[2]
            end
        end
    gmat = zeros(na,na,ny)
    trans = zeros(567,567)
    for j = 1:ny
       for k = 1:na
          gmat[k,ind_asg[k,j],j] = 1
       end
       trans[(j-1)*na+1:j*na,:] = kron(Piy[j,:]',gmat[:,:,j])
    end
    trans2=trans'
    probst = (1/(ny*na))*ones(ny*na,1)
    test=1
    while test > dist_tol
        probst1 = trans2*probst
        test = maximum(abs.(probst1-probst))
        probst = probst1
    end

    #   vectorize the decision rule to be conformable with probst
    #   calculate new aggregate capital stock  meanK


    kk=apol_egm[:]
    meank=probst'*kk

    # Using QuantEcon Invariant finder
    trans3 = trans
    for i = 1:size(trans,2)
        trans3[i,:] =trans[i,:]./sum(trans[i,:])
    end
    mc2 = MarkovChain(trans3)
    stationary_probs = stationary_distributions(mc2)[:, 1][1]
    meank_s = dot(stationary_probs,kk)

    # Calculate according to Muller
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
    A_supply2 = A_supply
    # GET the CDF
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
    ### SIMULATION
    T = 10^4
    Noind = 10^4
    at      = zeros(Noind,T+1)
    yt      = zeros(Noind,T)
    ct      = zeros(Noind,T)
    ht      = zeros(Noind,T)

    #at[:,1] = ones(Noind,1)           # initial asset level
    state_it = zeros(Noind,T)
    state_idx = zeros(Noind,T)
    t_s = [1/ny]
    for i =1:Noind
        #rng = MersenneTwister(1234)
        #s0 = rand(rng, 1)
        s0 = rand(1)
        s1 = (s0<=t_s)+(s0>t_s&&s0<=(t_s.*2)).*2+(s0>(t_s.*2)&&s0<=(t_s.*3)).*3+(s0>(t_s.*3)&&s0<=(t_s.*4)).*4+(s0>(t_s.*4)&&s0<=(t_s.*5)).*5+(s0>(t_s.*5)&&s0<=(t_s.*6)).*6+(s0>(t_s.*6)).*7
        state_it[i,:] = exp.(simulate(mc, T; init = s1))
        # Extract index
    end
    for i = 1:T
        #try to trace them nigrows
        #=
         ct[:,i] = (state_it[:,i].==egrid[1]).*c_interp(agrid,cpol_mat[:,1],at[:,i])+(state_it[:,i].==egrid[2]).*c_interp(agrid,cpol_mat[:,2],at[:,i])+(state_it[:,i].==egrid[3]).*c_interp(agrid,cpol_mat[:,3],at[:,i])                        +(state_it[:,i].==egrid[4]).*c_interp(agrid,cpol_mat[:,4],at[:,i])+(state_it[:,i].==egrid[5]).*c_interp(agrid,cpol_mat[:,5],at[:,i])+(state_it[:,i].==egrid[6]).*c_interp(agrid,cpol_mat[:,6],at[:,i])+(state_it[:,i].==egrid[7]).*c_interp(agrid,cpol_mat[:,7],at[:,i])
         # future assets
         at[:,i+1] = (1.0+r).*at[:,i]+state_it[:,i].*w-ct[:,i]
         #extract index
         =#
         ct[:,i] = (state_it[:,i].==egrid[1]).*c_interp(agrid,cpol_mat[:,1],at[:,i])+(state_it[:,i].==egrid[2]).*c_interp(agrid,cpol_mat[:,2],at[:,i])+(state_it[:,i].==egrid[3]).*c_interp(agrid,cpol_mat[:,3],at[:,i])                        +(state_it[:,i].==egrid[4]).*c_interp(agrid,cpol_mat[:,4],at[:,i])+(state_it[:,i].==egrid[5]).*c_interp(agrid,cpol_mat[:,5],at[:,i])+(state_it[:,i].==egrid[6]).*c_interp(agrid,cpol_mat[:,6],at[:,i])+(state_it[:,i].==egrid[7]).*c_interp(agrid,cpol_mat[:,7],at[:,i])
         # future assets
         at[:,i+1] = (1.0+r).*at[:,i]+state_it[:,i].*w.*H-ct[:,i]

    end
    K_simul = mean(mean(at[:,T-100:T]))
    return (at = at,mean_K= A_supply2, meank=meank, meank_s=meank_s, K_simul=K_simul, D_invariant=Λnm1_mat, pdf_a = pdf_a, cdf_a=cdf_a[1:end-1])
    #return (mean_K= A_supply2)
end

function c_interp(Amat,cpo,val)
    c_aux = interp(Amat, cpo)
    return c_aux.(val)
end
