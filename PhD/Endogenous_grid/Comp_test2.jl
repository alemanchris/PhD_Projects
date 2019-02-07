
# Create an Instance of a Household
am = Household(r = 0.03, w = 0.956)
# Do manually
r = 0.03
w = r_to_w(r)
@unpack a_vals, s_vals, u, Piy,z_vals,z_chain = am
setup_R!(am.R, a_vals, s_vals, r, w, u)

aiyagari_ddp = DiscreteDP(am.R, am.Q, am.β)

# Compute the optimal policy
results = solve(aiyagari_ddp, PFI)
# Obtain the values of the policy functions
@unpack z_size, a_size, n, a_min = am
a_star = reshape([a_vals[results.sigma[s_i]] for s_i in 1:n], a_size, z_size)
a_star_idx = reshape([results.sigma[s_i] for s_i in 1:n], a_size, z_size)
p_astar1 = polyfit(a_vals,a_star[:,1],1)
p_astar2 = polyfit(a_vals,a_star[:,2],1)
a_star1 = p_astar1(a_vals)
a_star2 = p_astar2(a_vals)
typeof(a_star[:,1])
#index_star1a = a_star1.<a_min
#index_star2a = a_star2.<a_min
#a_star1[index_star1a] = repeat([a_min],length(a_star1[index_star1a]),1)
#a_star2[index_star2a] = repeat([a_min],length(a_star1[index_star2a]),1)

a_STAR = [a_star1 a_star2]
c_star = zeros(a_size,z_size)
c_star2 = zeros(a_size,z_size)
#plot(a_vals, a_star)
for i = 1:z_size
    c_star[:,i] = (1+r)*a_vals.+ w*z_vals[i] .- a_star[:,i]
    c_star2[:,i] = (1+r)*a_vals.+ w*z_vals[i] .- a_STAR[:,i]
end
#plot(a_vals,c_star2)

#plot(a_vals,test1)
# Compute the stationary distribution:
####################################################33
# Sargets way
stationary_probs = stationary_distributions(results.mc)[:, 1][1]
K1 = dot(am.s_vals[:, 1], stationary_probs)
# Plot distribution
distrib = reshape(stationary_probs,a_size,z_size)
pdf_a = sum(distrib,dims=2)
plot(a_vals,pdf_a)
####################################################3
# Sargets way 2
gmat = zeros(a_size,a_size,z_size)
l_t = length(stationary_probs)
trans = zeros(l_t,l_t)
for j = 1:z_size
   for k = 1:a_size
      gmat[k,a_star_idx[k,j],j] = 1
   end
   trans[(j-1)*a_size+1:j*a_size,:] = kron(Piy[j,:]',gmat[:,:,j])
end
trans2=trans'
probst = (1/(z_size*a_size))*ones(z_size*a_size,1)
test=1
dist_tol = 1e-7
while test > dist_tol
    global probst1 = trans2*probst
    global test = maximum(abs.(probst1-probst))
    global probst = probst1
end

#   vectorize the decision rule to be conformable with probst
#   calculate new aggregate capital stock  meanK


kk=am.s_vals[:,1]
meank=probst'*kk
# plot the distributionn
distrib2 = reshape(probst,a_size,z_size)
pdf_b = sum(distrib2,dims=2)
plot(a_vals,pdf_b)
#
p_dert = polyfit(a_vals,pdf_b[:],6)
typeof(pdf_b[:])
fited = p_dert(a_vals)
plot(a_vals,fited)
# Calculate interpolations of the consumption policy rule
# Two Options:
# 1: Add 1 to the grid and get new interpolated policy rule
# 3: Just make a new finer grid and define the consumption policy interpolated on this on this grid
# Lets first Look at the discretized consumption policy rules (to many big jumps?)
##################################################################33
# Simulation
T = 10^4
Noind = 10^3
at      = zeros(Noind,T+1)
yt      = zeros(Noind,T)
ct      = zeros(Noind,T)
ct2      = zeros(Noind,T)
at2      = zeros(Noind,T+1)
#ht      = zeros(Noind,T)

at[:,1] = ones(Noind,1)           # initial asset level
state_it = zeros(Noind,T)
state_idx = zeros(Noind,T)
t_s = [1/z_size]
for i =1:Noind
    s0 = rand(1)
    s1 = (s0<=t_s)+(s0>t_s).*2
    state_it[i,:] = simulate(z_chain, T; init = s1)
end
#indicator = zeros(Noind,1)

for i = 1:T
     ct[:,i] = (state_it[:,i].==z_vals[1]).*c_interp(a_vals,c_star[:,1],at[:,i])+(state_it[:,i].==z_vals[2]).*c_interp(a_vals,c_star[:,2],at[:,i])
     at[:,i+1] = (1.0+r)*at[:,i].+state_it[:,i]*w.-ct[:,i]
     ct2[:,i] = (state_it[:,i].==z_vals[1]).*c_interp(a_vals,c_star2[:,1],at2[:,i])+(state_it[:,i].==z_vals[2]).*c_interp(a_vals,c_star2[:,2],at2[:,i])
     at2[:,i+1] = (1.0+r)*at2[:,i].+state_it[:,i]*w.-ct2[:,i]
     #indicator= at[:,i+1].>a_max_1
     #at[indicator,i+1].=a_max_1
end
K_simul = mean(mean(at[:,T-100:T]))
data = at[:,T-100:T]
histogram(data[:])


#######################################################333
# Piecewise Linear Approximation

a_ast=repeat(a_vals, 1, z_size) + c_star - repeat(z_vals*w, a_size, 1)/(1+r)
#a_ast2=repeat(a_vals, 1, z_size) + c_star2 - repeat(z_vals*w, a_size, 1)/(1+r)
#@unpack A_supply, iter, cdf_a, pdf_a, agrid_finer = compute_invariant(Household(),a_star)
@unpack A_supply3, iter, cdf_a, pdf_d  = compute_invariant(Household(),a_ast)
#@unpack A_supply, iter, cdf_a, pdf_a, agrid_finer = compute_invariant(Household(),a_ast2)

pdf_b =pdf_a[1:end-200]
help_1 = length(pdf_b)
agrid_finer2 = range(0,20,length=help_1)
plot(agrid_finer, cdf_a[1:end-1])
plot(agrid_finer, pdf_a)
norm_lis = pdf_b./(sum(pdf_b))
plot(agrid_finer2,norm_lis)
pdf_an = pdf_a./(sum(pdf_a))
plot(agrid_finer,pdf_an)


agrdi_finer3 = range(0,20,length=length(pdf_d))
plot(agrdi_finer3,pdf_d)

#################################
# My inverse

@unpack A_supply2, iter3, cdf_b, pdf_b = m_invariant(Household(),a_STAR)

@unpack na, agrid, agrid_finer= Household()
@unpack maxit, dist_tol, b = Household()
@unpack Pyinv, Piy, egrid, ny,a_min,a_max = Household()
Lambda0 = zeros(na*2,ny)
for j = 1:ny
    for k = 1:na*2
        Lambda0[k,j] = (agrid_finer[k]-a_min)/(a_max-a_min)*Pyinv[j]
    end
end
# Make it continous
#l0 = @(a,lambda) interp1(agrid2,lambda,a,'linear'); %  I can condition de densities
#l0(agrid_finer,lambda,a)
# update distribution for every pair
#crit1 = 10^(-8);

lambda = Lambda0
l1 = zeros(na*2,ny)
aprime2 = zeros(na*2,ny)
for i=1:ny
    for j=1:na*2
        #aux = interp(agrid,a_STAR[:,i],agrid_finer[j])
        #aux2 = interp1(kap,condecis(:,i),agrid_finer[j])
        aprox = i_STAR(agrid,a_STAR[:,i],agrid_finer[j])
        aprime2[j,i] = maximum([minimum([aprox,a_max]),a_min])
        #condecis2(j,i) = max(aux2,0);
    end
end
#maxiter = 10
#Inverse mapping
#invaprime(a,i) =interp(a_STAR[:,i],agrid_finer,a)
#a_star_int = zeros(length(agrid_finer),z_size)

int_star_1 = interp(agrid,a_STAR[:,1])
int_star_2 = interp(agrid,a_STAR[:,2])
a_star_int1 = int_star_1.(agrid_finer)
a_star_int2 = int_star_2.(agrid_finer)
a_star_int = [a_star_int1 a_star_int2]
iter3 = 0
dist1 = 1
l1 = zeros(na*2,ny)
#range  = [a_min,amax]
asol = invaprime(agrid_finer[100],a_star_int[:,2],agrid_finer)
while iter<2000 && dist1>1e-4
    global iter3 = iter3+1
        for k = 1:na*2
            for j = 1:ny
                for i = 1:ny
                    #%dert= invaprime(a,s,r,agrid2',condecis2,agrid2(k),i,wage),X0)
                    #asol = fzero(invaprime(a,s,r,agrid2',condecis2,agrid2(k),i,wage))
                    #asol = invaprime(agrid_finer[k],i)
                    asol = invaprime(agrid_finer[k],a_star_int[:,i],agrid_finer)
                    asol2= maximum([minimum([asol,a_max]),a_min])
                    l1[k,j] = l1[k,j]+Piy[i,j]*lo(agrid_finer,lambda[:,i],asol2)
                end
            end
        end


    global dist1 = maximum(maximum(abs.(lambda-l1)))
    #dist1 = norm(abs.(lambda-l1)
    global lambda = l1
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
cdf_b[length_aux+1] = cdf_b[length_aux]
#repeat the last value with the same so that the pdf of the last is zero
# GET THE PDF
pdf_b = diff(cdf_b)
plot(agrid_finer,pdf_b)
