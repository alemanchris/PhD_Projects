# Graphing

##

# Plots on the plane
using Plots
plotly() # if activated plots online on plotly
ts_length = 100
ϵ_values = randn(ts_length)
plot(ϵ_values, color="blue")
# Ploting histograms with Plots
histogram(ϵ_values)

using Plots
plot(rand(4,4))

using Plots
r_values=randn(1000)
histogram(r_values)
# Plots in the plane

using Gadfly
Gadfly.push_theme(:default)
plot([sin,cos],0,30)
plot(x=rand(20),y=rand(20),Geom.point)
plot(x=rand(20),y=rand(20),Geom.histogram)

# Ploting histograms
using Gadfly
e_values=randn(100)
e_values2=randn(100000)
plot(x=e_values2, Geom.histogram)
# Plot from dataframes,
plot(dataset("ggplot2", "diamonds"), x="Price", Geom.histogram)
####
