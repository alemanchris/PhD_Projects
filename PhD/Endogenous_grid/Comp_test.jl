
# Create an Instance of a Household
am = Household(a_max = 20.0, r = 0.03, w = 0.956)
# Do manually
r = 0.03
w = r_to_w(r)
@unpack a_vals, s_vals, u, Piy = am
setup_R!(am.R, a_vals, s_vals, r, w, u)

aiyagari_ddp = DiscreteDP(am.R, am.Q, am.β)

# Compute the optimal policy
results = solve(aiyagari_ddp, PFI)
# Obtain the values of the policy functions
@unpack z_size, a_size, n = am
a_star = reshape([a_vals[results.sigma[s_i]] for s_i in 1:n], a_size, z_size)
a_star_idx = reshape([results.sigma[s_i] for s_i in 1:n], a_size, z_size)
# Compute the stationary distribution:
# Sargets way
stationary_probs = stationary_distributions(results.mc)[:, 1][1]
K1 = dot(am.s_vals[:, 1], stationary_probs)

# Sargets way 2
gmat = zeros(a_size,a_size,z_size)
trans = zeros(400,400)
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
# Simulation
results.sigma
