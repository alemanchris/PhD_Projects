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


rand(1)
rand(2)
num=rand(1)
rand
