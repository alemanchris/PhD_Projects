compute_aiyagari(param())

ones(2,2)
compute_invariant(param())

@unpack mean_K,meank,meank_s,D_invariant, pdf_a, cdf_a = compute_invariant(param(r=0.036))
dert = [1 4;
        2 6;
        3 10]
0.96*(1+0.042)
(1-0.96)/(0.96)
wage*egrid[1]/0.02
alpha1 = 0.36
delta = 0.08
wage = (1-alpha1)*((alpha1/(0.02+delta))^alpha1)^(1/(1-alpha1))
rate = range(0.0, 0.02, length = 10)
mean_K_line = zeros(10,1)
for (i,v) in enumerate(rate)
    @unpack mean_K = compute_invariant(param(r=v))
    mean_K_line[i,1]=mean_K
end
plot(mean_K_line,rate)
@unpack cpol_mat, a_ast, apol_egm = compute_aiyagari(param(),plots=1)
@unpack cpol_mat, a_ast, apol_egm = compute_aiyagari(param(r= 0.036))

cpol_mat2 = compute_aiyagari(param())
compute_aiyagari(param(),plots=1)
@unpack agrid, agrid_finer,w = param(r=0.036)
plot(agrid, cpol_mat)
plot(agrid_finer, D_invariant)
plot(agrid_finer, pdf_a)
plot(agrid_finer, cdf_a)
@unpack egrid=rouwenhorst()

f(x)=compute_invariant(param(r=x))
function eq_r(x;) # pretending we don't know sqrt()
    mean_K = compute_invariant(param(r=x))
    z = mean_K
    return z
end
using Optim
using Roots
using Optim: converged, maximum, maximizer, minimizer, iterations #some extra functions
result = maximize(f, 0, 1)
converged(result) || error("Failed to converge in $(iterations(result)) iterations")
xmin = maximizer(result)
fmax = maximum(result)

fzero(eq_r, 0.04, 0.05)
dert2= f(0.04)
dert3= eq_r(0.05)
mean_K3=zeros(1,2)
dert = [1 2 3 4 5]'
dert2 = range(10,step=1,stop=14)
dot(dert2,dert)
sum(pdf_a)

dert = searchsortedfirst([1, 2, 4, 5, 14], 0)
typeof(dert)
