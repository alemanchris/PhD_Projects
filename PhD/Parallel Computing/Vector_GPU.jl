function solve_age_gpu(currentState::modelState_a)
    age     = currentState.age
    ne      = currentState.ne
    nx      = currentState.nx
    T       = currentState.T
    P       = currentState.P
    xgrid   = currentState.xgrid
    egrid   = currentState.egrid
    ssigma  = currentState.ssigma
    bbeta   = currentState.bbeta
    w       = currentState.w
    r       = currentState.r
    V       = currentState.V
    # Manipulate Matrices
    P_g = AFArray(permutedims(repeat(reshape(P',(ne,ne,1)), outer = [1,1,nx]),[3,2,1]))
    V_g = AFArray(permutedims(repeat(reshape(V,(nx,ne,1)), outer = [1,1,ne]),[1,3,2]))

    # Calculate expected utility
    expected = AFArray(zeros(nx,ne,nx))
    if (age < T)
        expected = AFArray(repeat(sum(P_g.*V_g,dims=3),outer = [1,1,nx]))
    end

    # Manipulate Matrices Budget constraint
    xgrid_g  = AFArray(permutedims(repeat(reshape(xgrid',(nx,1,1)),outer = [1,nx,ne]),[1,3,2]))
    xgridp_g = AFArray(permutedims(repeat(reshape(xgrid,(1,nx,1)),outer = [nx,1,ne]),[1,3,2]))
    egrid_g  = AFArray(repeat(reshape(egrid,(1,ne,1)),outer = [nx,1,nx]))

    # Compute consumption
    cons = (1+r).*xgrid_g +  egrid_g.*w - xgridp_g
    # Compute Utility
    utility = (cons.^(1-ssigma))/(1-ssigma) + bbeta.*expected
    # Eliminate negatives
    neg_ind = cons.<=AFArray(zeros(nx,ne,nx))
    utility = utility.*(.!neg_ind)+neg_ind.*(-10.0^(5))
    VV = Array(reshape(maximum(utility,dims = 3),(nx,ne))) #reshape is allowed

    return VV
end
function solve_VFI_gpu(T,ne,nx,P,xgrid,egrid,ssigma,bbeta,w,r,V_tomorrow_gpu,V_gpu)
    for age = T:-1:1
        currentState = modelState_a(ne,nx,age,T,P,xgrid,egrid,ssigma,bbeta,V_tomorrow_gpu,w,r)
        tempV_gpu = solve_age_gpu(currentState)
        V_gpu[age,:,:] = tempV_gpu
        V_tomorrow_gpu = tempV_gpu
    end

    return V_gpu
end
