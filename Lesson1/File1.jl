# Graphing

##

# Plots on the plane
using Plots
plotly() # if activated plots online on plotly
ts_length = 100
ϵ_values = randn(ts_length)
plot(ϵ_values, color="blue")


using Plots
plot(rand(4,4))
gui()
# Plots in the plane

using Gadfly
Gadfly.push_theme(:default)
plot([sin,cos],0,30)
plot(x=rand(20),y=rand(20),Geom.point)


####
x=10
y=linspace(2,100,30)
