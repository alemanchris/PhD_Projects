mcm = param()

    # retrieving parameters
@unpack beta, delta, na, agrid, agrid_finer, gamma, alpha1 = param()
@unpack r, w, maxit, dist_tol, b = param()

#@unpack Pyinv, Piy, egrid, ny, ygrid = rouwenhorst()
egrid = [0.301194211912202,0.449328964117222,0.670320046035639,1,1.49182469764127,2.22554092849247,3.32011692273655]
Pyinv = [0.00628217827163228,0.0608491084936224,0.241700981199883,0.382335464069726,0.241700981199883,0.0608491084936224,0.00628217827163229]
Piy =  [0.0262397497796233 0.152923483594817 0.361483063911416 0.328567584707170 0.114741751017859 0.0152660804766739 0.000778286512441828;
            0.0160443669891157 0.114741751017859 0.328567584707170 0.361483063911415 0.152923483594817 0.0247005562212874 0.00153919355833587;
            0.00945177150556068	0.0828345049143809 0.287445156270025 0.382789297365280 0.196113758787202 0.0384369616149092 0.00292854954264210;
            0.00536221880153007	0.0575309935179400 0.242023809517260 0.390165956326541 0.242023809517260 0.0575309935179400 0.00536221880153009;
            0.00292854954264206	0.0384369616149092 0.196113758787202 0.382789297365280 0.287445156270025 0.0828345049143808 0.00945177150556065;
            0.00153919355833585	0.0247005562212875 0.152923483594817 0.361483063911416 0.328567584707170 0.114741751017859 0.0160443669891157;
            0.000778286512441806 0.0152660804766738 0.114741751017859 0.328567584707170 0.361483063911416 0.152923483594817 0.0262397497796233]
#
ygrid = log.(egrid)
ny = 7
income_grid = egrid*w
@unpack cpol_mat, a_ast, apol_egm = compute_aiyagari(mcm)
    # Initial Values
Λ_invariant_init = zeros(length(agrid_finer), ny)
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


    #
# Discretization Sanrgent
# Getting the indices
ind_asg = Array{Int64}(undef,na,ny)
for i_z = 1:ny
    for i_a = 1:na
        ind_aux = minimum([searchsortedfirst(agrid, apol_egm[i_a,i_z]),na])
        ind_asg[i_a,i_z] = ind_aux
    end
end
#gmat=Array{Int64}(0,na,na,ny)
gmat = zeros(na,na,ny)
gmat[1,ind_asg[1,3],3]
trans = Array{Float64}(undef,567,567)
trans =zeros(567,567)
for j = 1:ny
    for k = 1:na
          gmat[k,ind_asg[k,j],j] = 1
    end
    trans[(j-1)*na+1:j*na,:] = kron(Piy[j,:]',gmat[:,:,j])
end

trans=trans'
probst = (1/(ny*na))*ones(ny*na,1)
test=1
while test > dist_tol
    global probst1 = trans*probst
    global test = maximum(abs.(probst1-probst))
    global probst = probst1
end
iteration=0
while iteration <10
    global iteration =iteration+1
end

kk = apol_egm[:]
meank=probst'*kk
#
A = [1 2;3 4]
B = [1 2;3 4]

kron(A,B)
dert = Array{Float64}(undef,5,5)
dert2 = repeat([1 2],4,5)
dert3 = zeros(2,2).+0.4
