mcm = param()

    # retrieving parameters
@unpack beta, delta, na, agrid, agrid_finer, gamma, alpha1 = param()
@unpack r, w, maxit, dist_tol, b = param()
@unpack Pyinv, Piy, egrid, ny, ygrid, mc = markov_aprox()

income_grid = egrid*w
@unpack cpol_mat, a_ast, apol_egm = compute_aiyagari(mcm,plots=1)

T = 3000
Noind = 3000
at      = zeros(Noind,T+1)
yt      = zeros(Noind,T)
ct      = zeros(Noind,T)
ht      = zeros(Noind,T)

at[:,1] = ones(Noind,1)           # initial asset level
state_it = zeros(Noind,T)

t_s = [1/ny]


#for t = 1:T
H = dot(egrid, Pyinv)
mpl_fun(K)  = (1.0-alpha1)*K^alpha1*H^-alpha1                        # computes w
inv_mpk(rate1) = H^(1.0-alpha1)*(alpha1/(rate1+delta))^(1.0/(1.0-alpha1))
KK0   = inv_mpk(r)
ww0   = mpl_fun(KK0)
for i =1:Noind
    s0 = rand(1)
    s1 = (s0<=t_s)+(s0>t_s&&s0<=(t_s.*2)).*2+(s0>(t_s.*2)&&s0<=(t_s.*3)).*3+(s0>(t_s.*3)&&s0<=(t_s.*4)).*4+(s0>(t_s.*4)&&s0<=(t_s.*5)).*5+(s0>(t_s.*5)&&s0<=(t_s.*6)).*6+(s0>(t_s.*6)).*7
    state_it[i,:] = exp.(simulate(mc, T; init = s1))
end
for i = 1:T
     ct[:,i] = (state_it[:,i].==egrid[1]).*c_interp(agrid,cpol_mat[:,1],at[:,i])+(state_it[:,i].==egrid[2]).*c_interp(agrid,cpol_mat[:,2],at[:,i])+(state_it[:,i].==egrid[3]).*c_interp(agrid,cpol_mat[:,3],at[:,i])                        +(state_it[:,i].==egrid[4]).*c_interp(agrid,cpol_mat[:,4],at[:,i])+(state_it[:,i].==egrid[5]).*c_interp(agrid,cpol_mat[:,5],at[:,i])+(state_it[:,i].==egrid[6]).*c_interp(agrid,cpol_mat[:,6],at[:,i])+(state_it[:,i].==egrid[7]).*c_interp(agrid,cpol_mat[:,7],at[:,i])
     # future assets
     at[:,i+1] = (1.0+r).*at[:,i]+state_it[:,i].*w-ct[:,i]
end
#=
for i = 1:T
     ct[:,i] = (state_it[:,i].==repeat([egrid[1]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,1],at[:,i])+(state_it[:,i].==repeat([egrid[2]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,2],at[:,i])+(state_it[:,i].==repeat([egrid[3]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,3],at[:,i])                        +(state_it[:,i].==repeat([egrid[4]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,4],at[:,i])+(state_it[:,i].==repeat([egrid[5]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,5],at[:,i])+(state_it[:,i].==repeat([egrid[6]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,6],at[:,i])+(state_it[:,i].==repeat([egrid[7]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,7],at[:,i])
     # future assets
     at[:,i+1] = (1+r).*at[:,i]+state_it[:,i].*w-ct[:,i]
end
=#
#=
for i = 1:T
     ct[:,t] = (state_it[:,i].==repeat([egrid[1]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,1],at[:,i])
                  +(state_it[:,i].==repeat([egrid[2]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,2],at[:,i])
                   +(state_it[:,i].==repeat([egrid[3]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,3],at[:,i])
                    +(state_it[:,i].==repeat([egrid[4]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,4],at[:,i])
                     +(state_it[:,i].==repeat([egrid[5]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,5],at[:,i])
                      +(state_it[:,i].==repeat([egrid[6]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,6],at[:,i])
                       +(state_it[:,i].==repeat([egrid[7]],inner=(Noind,1))).*c_interp(agrid,cpol_mat[:,7],at[:,i])
       # future assets
     at[:,i+1] = (1+r).*at[:,i]+state_it[:,i].*w-ct[:,i]
end
=#
K_simul = mean(mean(at[:,T-100:T]))
histogram(at[:,3001])
10^4
function c_interp(Amat,cpo,val)
    c_aux = LinearInterpolation(Amat, cpo)
    return c_aux(val)
end

repeat([1 2; 3 4], inner=(2, 1))
A = [1 2 3 4 5]
A.==2
dert = zeros(5,1)
rand(1)
dert[1,1]=rand(1)
for i =1:5
    #rng = MersenneTwister(1234)
    #s0 = rand(rng, 1)
    dert2 = rand(1)
    dert[i,1] = dert2

end
length(agrid)
plot(agrid,cpol_mat)


repeat(egrid', na, 1)*w*H
minimum(minimum(ct[:,T-100:T]))
