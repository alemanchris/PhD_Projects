# Excersice 1
function factorial2(n)
    comput = 1
    for i=1:n
        comput = comput*i
    end
    return comput
end
# Test factorial2
factorial(8)
factorial2(8)

# Excercise 2
function binomial_rv(n,p)
    suc = 0
    for i = 1:n
        num1 = rand(1)
        if num1[1]<p
            suc=suc+1
        end
    end
    return suc
end

binomial_rv(10000,0.5)
# Excercise 3 Montecarlo
count = 0
n = 1000000
for i =1:n
    # Estimate the are of a unit circle
    u = rand(2)
    if sqrt(u[1]^2 + u[2]^2)<=1
        count = count+1

    end
end
q_aprox = count/n #Aprox of B
pi_aprox = 4*q_aprox

# Excercise 4 Simulate and plot correlated time series

using Distributions
using Gadfly
alpha = 0.9
T     = 200
y = Array{Float64}(2)
x = Array{Float64}(T)
epsilon = randn(T)
x[1] = 0

for i = 2:T
    x[i] = alpha*x[i-1] + epsilon[i]
end

# Excercise 6 Plot the simulated series
using Gadfly
Gadfly.push_theme(:default)
alpha_2 = [0.0, 0.8, 0.98]
alpha_3 = [0.0 0.8 0.98] # see the difference with the one above
T = 200
y = Array{Float64}(T,3)
epsilon_2 = randn(T,3)
y[1,:] = 0
for j = 1:3
    for i = 2:T
        y[i,j] = alpha_2[j]*y[i-1,j] + epsilon_2[i,j]
    end
end
y1 = y[:,1]
y2 = y[:,2]
y3 = y[:,3]
x = 1:T
# default colors dont work on layers
plot(
  layer(x=x, y=y1, Geom.line, Theme(default_color=color("green"))),
  layer(x=x, y=y2, Geom.line, Theme(default_color=color("red"))),
  layer(x=x, y=y3, Geom.line, Theme(default_color=color("blue")))
)
