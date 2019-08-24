using InstantiateFromURL
activate_github("QuantEcon/QuantEconLecturePackages", tag = "v0.9.8");
using LinearAlgebra, Statistics, Compat
using Parameters, Plots, QuantEcon
using Interpolations
using StatsBase

# These functions are for Sargents QuantEcon
Household = @with_kw (r = 0.01,
                      w = 1.0,
                      σ = 1.0,
                      β = 0.96,
                      dist_tol = 1e-7,
                      # Three methods for discretization of the AR1
                      # 1. Markov chain: manually intorduce
                      z_chain = MarkovChain([0.9 0.1; 0.1 0.9], [0.1; 1.0]),
                      Pyinv = stationary_distributions(z_chain)[1],
                      # 2. Rouwenhorst:
                      #z_chain = rouwenhorst_mine(),
                      # 3. Tauchen:
                      #z_chain = markov_tauchen(),
                      a_min = 1e-10,
                      a_max = 18.0,
                      a_size = 200,
                      a_size_f = 300, #Finer Grid
                      a_vals   = range(a_min, a_max, length = a_size),
                      a_vals_f = range(a_min, a_max, length = a_size_f),
                      z_size = length(z_chain.state_values),
                      n = a_size * z_size,
                      s_vals = gridmake(a_vals, z_chain.state_values),
                      s_i_vals = gridmake(1:a_size, 1:z_size),
                      u = σ == 1 ? x -> log(x) : x -> (x^(1 - σ) - 1) / (1 - σ),
                      R = setup_R!(fill(-Inf, n, a_size), a_vals, s_vals, r, w, u),
                      # -Inf is the utility of dying (0 consumption)
                      Q = setup_Q!(zeros(n, a_size, n), s_i_vals, z_chain.p)) # Was just z_chain.p

struct pi_res_s
    # Generic structure to store results of Rouwenhorst or Tauchen
    Pyinv::Array{Float64}  # Invariant
    p::Array{Float64}       # Transition matrix
    egrid::Array{Float64} #egrid::Array{Float64}
    state_values::Array{Float64} #ygrid
    ny::Int64               # Number of states
    #pi_res(Pi_inv,Pi_y,egrid,ygrid,ny)
end

function rouwenhorst_mine(;
                    rho = 0.5,
                    p = (1.0+rho)/2, #p = (1.0+rho)/2,
                    q = p,
                    sigma_e = 0.4, # Dont know
                    sigma_y = sqrt((sigma_e^2)/(1-rho^2)),
                    ny = 7, # Number of states
                    my = 1) # mean of y

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
        pi_rouwenhorst = pi_res_s(Pyinv,Piy,egrid,ygrid,ny) #Save in structure form
        return pi_rouwenhorst
        #return (Pyinv=Pyinv, Piy=Piy ,egrid=egrid ,ny=ny, ygrid=ygrid)
        #return Pyinv
           # Persistence of the income process
end

function markov_tauchen(;
            ny = 7,
            rho = 0.2,
            sigmaint = 0.4,
            sigma    = sigmaint*sqrt(1-rho^2),
            mc = tauchen(ny, rho, sigma), # Generates same process as markovappr in Sargent
            Pyinv = stationary_distributions(mc)[1], #[1] is the method
            Piy = mc.p,
            ygrid = mc.state_values,
            egrid = exp.(ygrid))

    pi_tauchen = pi_res_s(Pyinv,Piy,egrid,ygrid,ny)
    return pi_tauchen
end

function setup_Q!(Q, s_i_vals, z_chain)
    for next_s_i in 1:size(Q, 3)
        for a_i in 1:size(Q, 2)
            for s_i in 1:size(Q, 1)
                z_i = s_i_vals[s_i, 2]
                next_z_i = s_i_vals[next_s_i, 2]
                next_a_i = s_i_vals[next_s_i, 1]
                if next_a_i == a_i
                    #Q[s_i, a_i, next_s_i] = z_chain.p[z_i, next_z_i]
                    Q[s_i, a_i, next_s_i] = z_chain[z_i, next_z_i]
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

# Firms' parameters
A = nothing
#const A = 1
A = 1
const N = 1
const α = 0.33
const β = 0.96
const δ = 0.05

function r_to_w(r)
    return A * (1 - α) * (A * α / (r + δ)) ^ (α / (1 - α))
end

function rd(K)
    return A * α * (N / K) ^ (1 - α) - δ
end

function solve_asset_policy_rule(am,r)
    # Set up problem
    w = r_to_w(r)
    @unpack a_vals, s_vals, u, a_size, z_size, n = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)

    aiyagari_ddp = DiscreteDP(am.R, am.Q, am.β)

    # Compute the optimal policy
    results = solve(aiyagari_ddp, PFI)

    # Get the asset policy function from results
    # Get all optimal actions across the set of a indices with z fixed in each column
    a_star = reshape([a_vals[results.sigma[s_i]] for s_i in 1:n], a_size, z_size)
    return a_star
end

function dist_sargent_A(am, r)

    # Set up problem
    w = r_to_w(r)
    @unpack a_vals, s_vals, u, a_size, z_size, n = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)

    aiyagari_ddp = DiscreteDP(am.R, am.Q, am.β)

    # Compute the optimal policy
    results = solve(aiyagari_ddp, PFI)

    # Get the asset policy function from results
    # Get all optimal actions across the set of a indices with z fixed in each column
    a_star = reshape([a_vals[results.sigma[s_i]] for s_i in 1:n], a_size, z_size)


    # Compute the stationary distribution
    stationary_probs = stationary_distributions(results.mc)[:, 1][1]

    # Calculate Aggregate Capital demand
    Kap = dot(am.s_vals[:, 1], stationary_probs)

    # Return K
    return (Kap_A = Kap, SSdist_A = stationary_probs)
end

# These are functions for manual discretization
function dist_sargent_B(am,r)
    # Set up problem (I should create a function)
    w = r_to_w(r)
    @unpack a_vals, s_vals, u, a_size, z_size, n, dist_tol, z_chain = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)

    a_star = solve_asset_policy_rule(am,r)

    na = a_size
    ny = z_size #2#z_chain.ny
    # Find the Indicees of a_star, on the grid, useful to find the "matrix" gmat
    ind_asg = Array{Int64}(undef,na,ny)
    for i_z = 1:ny
        for i_a = 1:na
            useles = findmin(abs.(a_vals.-a_star[i_a,i_z]))
              #useles = findmin(abs.(agrid.-apol_egm[i_a,i_z]))
              #ind_aux = minimum([searchsortedfirst(agrid, apol_egm[i_a,i_z]),na])
              ind_asg[i_a,i_z] = useles[2]
        end
    end
    gmat = zeros(na,na,ny)
    trans = zeros(na*ny,na*ny) # na*ny
    for j = 1:ny
        for k = 1:na
          gmat[k,ind_asg[k,j],j] = 1
        end
        #trans[(j-1)*na+1:j*na,:] = kron(Piy[j,:]',gmat[:,:,j])
        trans[(j-1)*na+1:j*na,:] = kron(am.z_chain.p[j,:]',gmat[:,:,j])
    end
    trans2 = trans'
    probst = (1/(ny*na))*ones(ny*na,1)
    test = 1
    while test > dist_tol
        probst1 = trans2*probst
        test = maximum(abs.(probst1-probst))
        probst = probst1
    end

    #   vectorize the decision rule to be conformable with probst
    #   calculate new aggregate capital stock  meanK


    #kk=apol_egm[:]
    kk = a_star[:]
    Kap=probst'*kk

    # Using QuantEcon Invariant finder
    #=
    trans3 = trans
    for i = 1:size(trans,2)
        trans3[i,:] =trans[i,:]./sum(trans[i,:])
    end
    mc2 = MarkovChain(trans3)
    stationary_probs = stationary_distributions(mc2)[:, 1][1]
    Kap = dot(stationary_probs,kk) #must be wrong
    =#
    return (Kap_B = Kap, SSdist_B = probst)
end

# 3. Simulation

function simulation_C(am,r)
    # Set up problem (I should create a function)
    w = r_to_w(r)
    @unpack a_vals, s_vals, u, a_size, z_size, n, dist_tol, z_chain = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)
    na = a_size
    ny = z_size #2#z_chain.ny

    a_star = solve_asset_policy_rule(am,r)
    c_star = zeros(a_size,z_size)

    for i = 1:ny
        c_star[:,i] = (1+r)*a_vals.+ w*z_chain.state_values[i] .- a_star[:,i]
    end
    T = 10^4 # if 5 then out of memory
    Noind = 10^4
    at      = zeros(Noind,T+1)
    yt      = zeros(Noind,T)
    ct      = zeros(Noind,T)
    #ht      = zeros(Noind,T)

    #at[:,1] = ones(Noind,1)           # initial asset level
    state_it = zeros(Noind,T)
    state_idx = zeros(Noind,T)
    t_s = [1/ny]
    for i =1:Noind
        #rng = MersenneTwister(1234)
        #s0 = rand(rng, 1)
        s0 = rand(1)
        #s1 = (s0<=t_s)+(s0>t_s&&s0<=(t_s.*2)).*2+(s0>(t_s.*2)&&s0<=(t_s.*3)).*3+(s0>(t_s.*3)&&s0<=(t_s.*4)).*4+(s0>(t_s.*4)&&s0<=(t_s.*5)).*5+(s0>(t_s.*5)&&s0<=(t_s.*6)).*6+(s0>(t_s.*6)).*7
        s1 = (s0<=t_s)+(s0>(t_s)).*2
        state_it[i,:] = simulate(z_chain, T; init = s1)
    end
    for i = 1:T
        #try to trace them nigrows
        #
         ct[:,i] = (state_it[:,i].==z_chain.state_values[1]).*val_interp(a_vals,c_star[:,1],at[:,i])+(state_it[:,i].==z_chain.state_values[2]).*val_interp(a_vals,c_star[:,2],at[:,i])
         # future assets
         at[:,i+1] = (1.0+r).*at[:,i]+state_it[:,i].*w-ct[:,i]
         #extract index
         #
        #= Variant with no H (But with more income states, find a way to loopit)
         ct[:,i] = (state_it[:,i].==egrid[1]).*c_interp(agrid,cpol_mat[:,1],at[:,i])+(state_it[:,i].==egrid[2]).*c_interp(agrid,cpol_mat[:,2],at[:,i])+(state_it[:,i].==egrid[3]).*c_interp(agrid,cpol_mat[:,3],at[:,i])                        +(state_it[:,i].==egrid[4]).*c_interp(agrid,cpol_mat[:,4],at[:,i])+(state_it[:,i].==egrid[5]).*c_interp(agrid,cpol_mat[:,5],at[:,i])+(state_it[:,i].==egrid[6]).*c_interp(agrid,cpol_mat[:,6],at[:,i])+(state_it[:,i].==egrid[7]).*c_interp(agrid,cpol_mat[:,7],at[:,i])
         # future assets
         at[:,i+1] = (1.0+r).*at[:,i]+state_it[:,i].*w-ct[:,i]
         #extract index
         =#
        #= Variant with H (incomplete)
         ct[:,i] = (state_it[:,i].==egrid[1]).*c_interp(agrid,cpol_mat[:,1],at[:,i])+(state_it[:,i].==egrid[2]).*c_interp(agrid,cpol_mat[:,2],at[:,i])+(state_it[:,i].==egrid[3]).*c_interp(agrid,cpol_mat[:,3],at[:,i])                        +(state_it[:,i].==egrid[4]).*c_interp(agrid,cpol_mat[:,4],at[:,i])+(state_it[:,i].==egrid[5]).*c_interp(agrid,cpol_mat[:,5],at[:,i])+(state_it[:,i].==egrid[6]).*c_interp(agrid,cpol_mat[:,6],at[:,i])+(state_it[:,i].==egrid[7]).*c_interp(agrid,cpol_mat[:,7],at[:,i])
         # future assets
         at[:,i+1] = (1.0+r).*at[:,i]+state_it[:,i].*w.*H-ct[:,i]
         =#
    end
    K_simul = mean(mean(at[:,T-100:T]))
    # Obtain weights
    data_vec = vec(at[:,T-100:T]) #Vectorize data
    bin_it = fit(Histogram,data_vec; closed=:left,nbins=na)
    return (Kap_C = K_simul, SSdist_C = bin_it.weights)
#return (mean_K= A_supply2)
end

function val_interp(Amat,cpo,val)
    # Returns interpolated value C at val
    c_aux = interp(Amat, cpo)
    return c_aux.(val)
end

# 4. Piece Wise Linear Interpolation (Failed attempt: Think again to find mistake)
function PWLinear_D(am,r)
    # Set up problem (I should create a function)
    w = r_to_w(r)
    @unpack a_vals, s_vals, u, a_size, z_size, n, dist_tol, z_chain, Pyinv, a_vals_f, a_size_f = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)
    na = a_size
    ny = z_size #2#z_chain.ny

    a_star = solve_asset_policy_rule(am,r)
    # find the associated endogenous grid, the value of assets today
    a_ast_mat = zeros(a_size,z_size)
    for i =1:ny
        a_ast_mat[:,i] = val_interp(a_star[:,i],a_vals,a_vals) # Todays Assets
    end
    # Initial guess one
#    Linit=ones(a_size_f,z_size)./sum(sum(ones(a_size_f,z_size))
    # Initial guess two (this is just a line)
    Linit = zeros(a_size_f,z_size)
    for i_y = 1:ny
        for (i_a, a_v) in enumerate(a_vals_f)
            Linit[i_a, i_y] = (a_v - a_vals[1]) / (a_vals[end] - a_vals[1]) * Pyinv[i_y]
        end
    end
    #
    Λn_mat   = zeros(a_size_f, ny)
    Λnm1_mat = copy(Linit)

    iter = 0
    dist = 10.0


    a_ast_itp = LinearInterpolation((a_vals, z_chain.state_values), a_ast_mat) # next period's assets, current y
    # find the associated endogenous grid, that is, value of assets today
    while iter<20000 && dist>1e-5
        iter += 1
        for i_y = 1:ny # next period
            for (i_a, a_v) in enumerate(a_vals_f) # next period
                vec_temp = zeros(ny)
                for (i_y0, y0_v) in enumerate(z_chain.state_values) # last period y
                    aux_ast = a_ast_itp(a_v, y0_v)
                    aval  = minimum([maximum([aux_ast, a_vals[1]]), a_vals[end]]) # today's assets (endogenous grid)
                    aux_ast_pos = searchsortedfirst(a_vals_f, aval)
                    #length_aux = length(ay.agrid_finer)
                    ind_r = minimum([maximum([aux_ast_pos, 2]), a_size_f]) # No [1] cuz if not, in following I will have ind_r-1 = 0
                    Λval = Λnm1_mat[ind_r-1, i_y0] + (Λnm1_mat[ind_r, i_y0]- Λnm1_mat[ind_r-1, i_y0]) / (a_vals_f[ind_r] - a_vals_f[ind_r-1]) * (aval - a_vals_f[ind_r-1])
                    vec_temp[i_y0] = Λval
                end
                Λn_mat[i_a, i_y] = dot(z_chain.p[:, i_y], vec_temp)
                #Λn_mat[i_a, i_y] = dot(z_chain.p[i_y, :]', vec_temp)
            end
        end
        dist = maximum(abs.(Λn_mat - Λnm1_mat))
        #copyto!(Λnm1_mat, Λn_mat)
        Λn_mat = copy(Λnm1_mat) # This is the CDF
        #Λnm1_mat .= Λn_mat
    end
    A_supply = 0.0

    for i_y = 1:ny

        sum_val_a = 0.0

        for i_a = 1:(a_size_f-1)
            anp1 = a_vals_f[i_a+1]
            an   = a_vals_f[i_a]
            sum_val_a += 0.5*(Λn_mat[i_a+1, i_y] - Λn_mat[i_a, i_y]) * (anp1 + an)
        end

        A_supply += sum_val_a + Λn_mat[1, i_y]*a_vals_f[1]
    end
    Kap_D = A_supply
    CDF_lambda_inv = Λn_mat
    CDF_inv = sum(CDF_lambda_inv,dims=2)
    PDF_inv = diff(CDF_inv,dims=1)
    return (Kap_D = Kap_D ,SSdist_D =PDF_inv)
end

# Piece Wise linear Interpolation of CDF (Success!)
function CDF_pwise_E(am,r)
    # Set up problem (I should create a function)
    w = r_to_w(r)
    @unpack a_vals, s_vals, u, a_size, z_size, n, dist_tol, z_chain, Pyinv, a_vals_f, a_size_f = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)
    na = a_size
    ny = z_size #2#z_chain.ny
    nk = a_size_f
    a_star = solve_asset_policy_rule(am,r)

    # Initialize distribution function
    nk0 = sum(a_vals_f.<=0)
    gk  = zeros(nk,ny)
    gk[nk0+1:nk,:] = ones(nk-nk0,ny)
    gk = gk.*Pyinv'
    # Compute the inverse
    # Elimination of the rows with the same entry as the subsequent row
    concat1 = [a_star[2:na,1];0]
    concat2 = [a_star[2:na,2];0]
    log_del_1 = a_star[:,1].==concat1
    log_del_2 = a_star[:,2].==concat2
    # Start eliminating (VERY ADHOCK TILL NOW)
    aopt_aux11 = a_vals[log_del_1]
    aopt_aux12 = a_star[log_del_1,1]
    aopt_aux21 = a_vals[log_del_2]
    aopt_aux22 = a_star[log_del_2,2]
    #=
    aopt_aux11 = copy(a_vals)
    aopt_aux12 = copy(a_star[:,1])
    aopt_aux21 = copy(a_vals)
    aopt_aux22 = copy(a_star[:,2])
    deleteat!(vec(aopt_aux11),vec(log_del_1))
    deleteat!(vec(aopt_aux12),vec(log_del_1))
    deleteat!(vec(aopt_aux21),vec(log_del_2))
    deleteat!(vec(aopt_aux22),vec(log_del_2))
    =#
    aopte = [aopt_aux11 aopt_aux12]
    aoptu = [aopt_aux21 aopt_aux22]
    a1 = zeros(nk,ny)
    for l = 1:ny
        for i = 1:nk
            if l==1
                if a_vals_f[i] <= aopte[1,2]  # a'> ag[i] for all a
                    a1[i,l] = a_vals[1]-0.1
                elseif a_vals_f[i]>=maximum(aopte[:,2])  #   ag[i]>a' for all a
                    a1[i,l] = maximum(aopte[:,1])
                else
                    a1[i,l] = val_interp(aopte[:,2],aopte[:,1],a_vals_f[i])
                end
            else
                if a_vals_f[i] <= aoptu[1,2]
                    a1[i,l] = a_vals[1]-0.1
                elseif a_vals_f[i]>=maximum(aoptu[:,2])
                    a1[i,l]=maximum(aoptu[:,1])
                else
                    a1[i,l]=val_interp(aoptu[:,2],aoptu[:,1],a_vals_f[i])
                end
            end
        end
    end
   # computation of invariant distribution of wealth..";
   q1 = 0
   ngk = 25000
   while  q1<ngk #until (q1>ngk);
       q1 = q1+1

       gk0 = copy(gk)
       gk = zeros(nk,ny)

       for l=1:ny
           for i=1:nk
               k0 = a1[i,l]
               for l1 = 1:ny
                   if k0<=a_vals_f[1]
                       gk[i,l1] = gk[i,l1]+0
                   elseif k0>=a_vals_f[nk]
                       gk[i,l1] = gk[i,l1]+z_chain.p[l,l1]*Pyinv[l]
                   else
                       gk[i,l1] = gk[i,l1]+val_interp(a_vals_f,gk0[:,l],k0)*z_chain.p[l,l1]
                   end
               end
           end
       end


       gk = gk./(sum(gk[nk,:]))
       kritg = maximum(abs.(gk0-gk))
   end
   # computing the mean
   kk1 = (gk[1,1]+gk[1,2])*a_vals_f[1]
   gk1 = (gk[2:nk,1]+gk[2:nk,2])-(gk[1:nk-1,1]+gk[1:nk-1,2])
   ag1 = (a_vals_f[2:nk]+a_vals_f[1:nk-1])./2
   kk1 = kk1 + gk1'*ag1

   ag1 = [a_vals_f[1];ag1]
   gk1 = [gk[1,:];(gk[2:nk,:]-gk[1:nk-1,:])]
   CDF_inv = sum(gk1,dims=2)
   PDF_inv = diff(CDF_inv,dims=1)
   KapE = kk1

    return (Kap_E = kk1, SSdist_E = PDF_inv)
end
