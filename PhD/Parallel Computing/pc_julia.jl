using Distributed
addprocs(4)
# Serial, Sequencial
function Value(age, ix, ie)
    VV = -10^3;
    for(ixp = 1:nx)
        expected = 0.0;
        if(age < T)
            for(iep = 1:ne)
                expected = expected + P[ie, iep]*V[age+1, ixp, iep];
            end
        end
        cons = (1 + r)*xgrid[ix] + egrid[ie]*w - xgrid[ixp];
        utility = (cons^(1-ssigma))/(1-ssigma) + bbeta*expected;
        if(cons <= 0)
            utility = -10^5;
        end
        if(utility >= VV)
            VV = utility;
        end
    end
    return(VV);
end

for(age = T:-1:1)
    for(ix = 1:nx)
        for(ie = 1:ne)
            V[age, ix, ie] = Value(age, ix, ie);
        end
    end
end
