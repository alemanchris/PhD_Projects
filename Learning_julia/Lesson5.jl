# Law of Large Numbers (LLN) and Central Limit theorem (CLT)
# Excersise 1
using InstantiateFromURL
# load common packages# activate the QuantEcon environment
activate_github("QuantEcon/QuantEconLecturePackages", tag = "v0.9.4");

using LinearAlgebra, Statistics, Compat
using StatPlots, Distributions, Random, Statistics

#using Plots, Distributions, Random, Statistics
plotly()
#gr(fmt = :png, size = (900, 500))
function exercise1(distribution = Uniform(0, π/2); n = 250, k = 10_000, g = sin, g′ = cos)
    μ, σ = mean(distribution), std(distribution)
    y = rand(distribution, n, k)
    y = mean(y, dims = 1)
    y = vec(y)
    error_obs = sqrt(n) .* (g.(y) .- g.(μ))
    density(error_obs, label = "Empirical Density")
    return plot!(Normal(0, g′(μ) .* σ), linestyle = :dash, label = "Asymptotic",
                 color = :black)
    #return gp
    #return plot(x=Normal(0, g′(μ) .* σ), Geom.line)
    #return μ, σ
end

exercise1()
