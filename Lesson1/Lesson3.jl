# Exercise 1
A = [0.8 -0.2;-0.1 0.7]
sigma = [0.5 0.4;0.4 0.6]
S0 = zeros(2,2)
function compute_wea(A,sigma,S0)
	Q = sigma*sigma'
	S = S0
	for i in 1:500
		S = A*S*A' + Q
	end
	return S
end

value = compute_wea(A,sigma,S0)
println(value)
using QuantEcon

solving = solve_discrete_lyapunov(A,sigma*sigma')
println(solving)
