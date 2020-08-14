##
function test()
    sum = 1.0
    for n = 1:20
        sum += (-7)^n / factorial(n)
    end
    return sum
end

test()

##
exp(-7)

##
function test()
    
    sum = 1.0
    for n = 1:20
        sum += 7^n / factorial(n)
    end
    return 1 / sum
end

test()

##
sum = 1.0
for n = 1:20
    global sum
    sum += 7^n / factorial(n)
end

1/ sum

##
n = 1
while n < 10
    global n
    n += 1
end