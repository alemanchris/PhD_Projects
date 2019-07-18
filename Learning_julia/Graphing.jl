# installing IPython
watch tutorial
# installing IJulia
check this folder and delete the images in folder julia-1.0
C:\Users\Aleman\AppData\Roaming\jupyter\kernels\julia-1.0
# follow this steps
# set the enviroment where you have installed the jupyter
ENV["JUPYTER"] = "C:\\Users\\Aleman\\AppData\\Local\\Programs\\Python\\Python37\\Scripts\\jupyter.exe"
# then
] add IJulia
# then 
] build IJulia
using IJulia
notebook()

if there is IJulia already the type
] rm IJulia

# Graphing
# To clear all Vars
Ctr+D to restart Julia
Press enter to star again
clearconsole()
# to open console
ctrl+j then Ctrl+o
to clear console
ctrl+j then ctrl+c
Open repl
Ctrl+j Ctrl+r - Open a REPL
Inline Evaluation
You can evaluate your Julia code inline by navigating your cursor to the appropriate code and hitting Ctrl+Enter. This will run the code block that the cursor is contained in. For example, if you go to the top of a for loop, it will evaluate the whole for loop, or if it's inside of a function, it will evaluate the function (i.e. create the function). To specifically choose which code to evaluate, highlight the appropriate parts and use Ctrl+Enter. To evaluate the whole script, use Ctrl+Shift+Enter.

##

# Plots on the plane
using Plots
plotly() # if activated plots online on plotly
ts_length = 100
ϵ_values = randn(ts_length)
plot(ϵ_values, color="blue")

# Ploting histograms with Plots
using Plots
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
plot(x=rand(20),y=rand(20),Geom.point) #works if run from console
plot(x=rand(20),y=rand(20),Geom.histogram)

# Ploting histograms
using Gadfly
Gadfly.push_theme(:default)
e_values=randn(100)
e_values2=randn(100000)
plot(x=e_values2, Geom.histogram)
# Plot from dataframes,
plot(dataset("ggplot2", "diamonds"), x="Price", Geom.histogram)
####

# Ploting multiple lines with layers
plot(layer( x=[1:10], y=rand(10),Geom.point, Geom.line, Theme(default_color=color("orange")) ),
      layer( x=[1:10], y=rand(10),Geom.point, Geom.line, Theme(default_color=color("purple"))) )
