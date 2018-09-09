# Text arrays and touples
actions = ["surf","ski"]
actions_2 = ( "surf" , "ski")
# some logicals not to confuse
#=

P && Q is true if both are true, otherwise it’s false
P || Q is false if both are false, otherwise it’s true
=#
# Exercise 1
# Part 1
# the use of zip()
x_vals = [1,2,3]
y_vals = [1,1,1]
sum([ i[1]*i[2] for i in zip(x_vals,y_vals)])
# or
sum([x*y for (x,y) in zip(x_vals, y_vals)])
# Part 2
# The use of a comprehension
count_s = 0
sum([ i%2==0 for i in 0:99 ])

# with out using comprehension
doubles = [2i for i in 1:4]
for i = 0:2:t
    y = y+1
end
# Part 3
pairs = [(2, 5), (4, 2), (9, 8), (12, 10)]
sum([x%2==0 &&  y%2==0 for (x,y) in pairs])

# Excersice 2
countries = ("Japan", "Korea", "China")
cities = ("Tokyo", "Seoul", "Beijing")
for (i, country) in enumerate(countries)
    city = cities[i]
    println("The capital of $country is $city")
end
function p(x,coef)
    value = 0
    for (i,a_s) in enumerate(coef)
        value = value+(a_s*x^(i-1))
    end
    return value
end

p(2,[1,2,3])
p(1,[2,4])
# The short version
p_2(x, coeff) = sum([a * x^(i-1) for (i, a) in enumerate(coeff)])
# Excersice 3
function num_uper(x)
    count = 0
    for i in x
        if i == uppercase(i) && isalpha(i)
            count = count+1
        end
    end
    return count
end
 num_uper("U RaN iG ro")
 # one liner
count_num(x) = sum([ j == uppercase(j) && isalpha(j) for j in x ])
count_num("DerTderTCD")
count_num("U RaN iG ro")
# Excersice 4
f_ex4_2(seq_a, seq_b) = issubset(Set(seq_a), Set(seq_b))

# Excersice 5
function linapprox(f, a, b, n, x)
    #=
    Evaluates the piecewise linear interpolant of f at x on the interval
    [a, b], with n evenly spaced grid points.

    =#
    length_of_interval = b - a
    num_subintervals = n - 1
    step = length_of_interval / num_subintervals

    # === find first grid point larger than x === #
    point = a
    while point <= x
        point += step
    end

    # === x must lie between the gridpoints (point - step) and point === #
    u, v = point - step, point

    return f(u) + (x - u) * (f(v) - f(u)) / (v - u)
end
f_ex5(x) = x^2
g_ex5(x) = linapprox(f_ex5, -1, 1, 3, x)
x_grid = linspace(-1, 1, 100)
y_vals = map(f_ex5, x_grid)
y_approx = map(g_ex5, x_grid)
using Gadfly
Gadfly.push_theme(:default)
plot(layer( x=x_grid, y=y_vals, Geom.line, Theme(default_color=color("orange")) ),
     layer( x=x_grid, y=y_approx, Geom.line, Theme(default_color=color("purple"))) )

# Excersice 6

f = open("us_cities.txt","r")
total_pop = 0
for line = eachline(f)
    city, population = split(line,':')
    total_pop += parse(Int, population)
end

f_ex6 = open("us_cities.txt", "r")
total_pop = 0
for line in eachline(f_ex6)
    city, population = split(line, ':')            # Tuple unpacking
    total_pop += parse(Int, population)
end
