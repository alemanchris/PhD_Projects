using Plots
plot(rand(4,4))


r_values=randn(1000)
histogram(r_values)

] Pkg.dir("WebIO", "assets")
;jupyter labextension install webio
;jupyter labextension enable webio/jupyterlab_entry

using Interact
ui = button()
display(ui)


import RDatasets
iris = RDatasets.dataset("datasets", "iris")
using StatPlots, Interact
using Blink
w = Window()
body!(w, dataviewer(iris))
