function solve_age_y(currentState::modelState_yy)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P       = currentState.P
    ut_g   = currentState.ut_b
    neg_ind_g   = currentState.neg_ind_b
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    #P_g = permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1])
    #V_g = permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2])
    for ie = 1:ne
        # Calculate expected utility
        #expected = zeros(nx,ne,nx)
        expected = zeros(nx,1)
        if (age < T)
            #expected = repeat(sum(P[ie,:].*V,dims=2),outer = [1,1,nx])
            expected = reshape(sum(P[ie,:].*V,dims=2),(nx,1))
        end

        # Manipulate Matrices Budget constraint
        #xgrid_g  = permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2])
        #xgridp_g = permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2])
        #egrid_g  = repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx])

        # Compute consumption
        #cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
        # Compute Utility
        utility = reshape(ut_g[:,ie,:],(nx,nx)) .+ bbeta.*expected
        # Eliminate negatives
        #neg_ind = cons.<=zeros(nx,ne,nx)
        utility = utility.*(.!reshape(neg_ind_g[:,ie,:],(nx,nx))).+reshape(neg_ind_g[:,ie,:],(nx,nx)).*(-10.0^(5))
        VV[:,ie] = reshape(maximum(utility,dims = 2),(nx,1)) #reshape is allowed
    end
    return VV
end

# Some funny test

dert1 = AFArray([1 2 3;4 5 6;7 8 9;10 11 12])
dert2 = [10 20 30]
dert3 = [10;20;30]

dert1.*dert1[1,:]
dert1.*reshape(dert1[1,:],(1,3))
dert1[1,:]

nigow = maximum(dert1[1,:])
typeof(nigow)

ind = nx*ne
ne = 15
ix      = convert(Int, floor((ind-0.05)/ne))+1
ie      = convert(Int, floor(mod(ind-0.05, ne))+1)
ie = convert(Int, floor((ind-0.05)/nx))+1
ix      = convert(Int, floor(mod(ind-0.05, nx))+1)
