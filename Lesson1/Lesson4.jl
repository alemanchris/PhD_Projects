using Distributions
mutable struct AR1
    a
    b
    σ
    ϕ
end
m = AR1(0.9, 1, 1, Beta(5, 5))
function simulate(m::AR1, n::Integer, x0::Real)
    x_2 = Array{Float64}(n)
    x_2[1] = x0
    for i in 2:n
        x_2[i] = m.a*x_2[i-1] + m.b + m.σ*rand(m.ϕ)
    end
    return x_2
end
series = simulate(m,100,3)
using Gadfly
Gadfly.push_theme(:default)
x = linspace(1,100,100)
plot(x=x, y=series, Geom.line)
typeof(series)

dert = rand(100)
