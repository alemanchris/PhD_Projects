# Julia Lesson 3 summary
using InstantiateFromURL
# activate the QuantEcon enviroment
activate_github("QuantEcon/QuantEconLecturePackages", tag = "v0.9.5")
# Load common packages
using LinearAlgebra, Statistics, Compat
a = [1,2,3] # comma is like semicolom
# Types look like this Array{Float64,1} meaning one dimension
typeof(randn(100))
# for the number of dymensions
ndims(a)
# to see the number of elements
size(a)
# A vector is a one dimensional array
# A matrix is a two dimensional array
# vectors are column arrays
# row arrays are matrices
# fill is like the repmat
fill(a,2,2)
# Create empty arrays
x = Array{Float64}(undef,2,2)
# Copy
x = [1,2,3]
y = copy(x)
# Similar is just the same size
z = similar(x)
# randn()
# Views and Slices
a = [1 2; 3 4]
b = a[:,2] # : is like copy. SHOW just shows the resut
# view does not copy the value
@views b = a[:,2]
@show b
a[:,2] = [4,5]
@show a
@show b
#
