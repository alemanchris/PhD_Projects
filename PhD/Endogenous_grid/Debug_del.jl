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
#=
    while  q1<ngk #until (q1>ngk);
       global q1 = q1+1
       global gk
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
=#
# computing the mean
kk1 = (gk[1,1]+gk[1,2])*a_vals_f[1]
gk1 = (gk[2:nk,1]+gk[2:nk,2])-(gk[1:nk-1,1]+gk[1:nk-1,2])
ag1 = (a_vals_f[2:nk]+a_vals_f[1:nk-1])./2
kk1 = kk1 + gk1'*ag1

ag1 = [a_vals_f[1];ag1]
french = (gk[2:nk,:]-gk[1:nk-1,:])
french3 = gk[1,:]'
gk1 = [gk[1,:]';(gk[2:nk,:]-gk[1:nk-1,:])]
CDF_inv = sum(gk1,dims=2)
PDF_inv = diff(CDF_inv,dims=1)
KapE = kk1

DERT = a_vals_f'

amin10=-2                 #/* asset grid */
amax10=3000
na0=201
astep0=(amax10-amin10)/(na0-1)
#a0=seqa(amin10,astep0,na0);
nk0=3*na0 #            /* asset grid for distribution */
agstep3=(amax10-amin10)/(nk0-1)

#
A_test = [1 2;3 4;5 6;7 8]
ldel1 = [2;2;2;2]
ldel2 = [0;2;0;2]
log_stuf = .!(ldel1.==ldel2)
new_dert2 = [A_test[log_stuf,1] A_test[log_stuf,2]]
inv_log = log_stuf'
new_dert3 = [A_test[inv_log,1] A_test[inv_log,2]]
something_again = A_test[inv_log,1]
