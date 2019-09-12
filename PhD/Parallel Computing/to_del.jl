

# Manipulate Matrices
P_g = permutedims(repeat(P',outer = [1,1,nx]),[3,2,1])
V_g = permutedims(repeat(V,outer = [1,1,ne]),[1,3,2])

# Calculate expected utility
expected = zeros(nx,ne,nx)
if (age < T)
    expected = repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx])
end

# Manipulate Matrices Budget constraint
xgrid_g  = permutedims(repeat(xgrid',outer = [1,nx,ne]),[1,3,2])
xgridp_g = permutedims(repeat(xgrid,outer = [nx,1,ne]),[1,3,2])
egrid_g  = repeat(egrid,outer = [nx,1,nx])

# Compute consumption
cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
# Compute Utility
utility = (cons.^(1-ssigma))/(1-ssigma) + bbeta.*expected
# Eliminate negatives
neg_ind = cons.<=zeros(nx,ne,nx)
utility[neg_ind] = neg_ind.*(-10.0^(5))
