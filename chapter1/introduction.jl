2 + 3

sin(pi / 4)

log(2, 4)

# column vector/ Int64
x = [10, 20, 30]

# Float 64
x = [20, 30, 40, 0.1]

# row vector
x = [20 30 40]

# transpose
x'

x = [10 * i for i = 1:5]

last(x)

minimum(x)

sum(x)

# ! means inplace?
append!(x, 99)

x

length(x)

# broadcasting
sin.(x)


##
using PyPlot
x = range(0, stop=2 * pi, length=100)
y = sin.(3x);
plot(x, y, color="red")
title("The sine function")
gcf() # this is needed to show plots https://github.com/julia-vscode/julia-vscode/issues/325

# matrix
A = [-1 0.2 0.7; 0.09 -1 0.2; 1 1 1]

A'

inv(A)

# almost identity
A * inv(A)

v = [0 0 1]

# this does not work
# A * v

A * v'

# solve linear eq
A \ v'

# identical
inv(A) * v'

A^2

A * A

# todo: syntax for element-wise power etc.

2 == 3

2 <= 3

(2 == 2) || (1 < 0)

(2 == 2) && (1 < 0)

iseven(2)

isodd(3)

# defining functions

function square_(x)
    return x^2
end

square_(3)

square_.([1 2])

cube_(x) = x^3

cube_(5)

# anonymous function
x -> x^3

filter(x -> x > 0, [-2, 3, 1, -1])

count(x -> x > 0, [-1, 3, 1, -1])

# typed function

function typed_square_(x::Float64)
    return x^2
end

typed_square_(5.5)

# typed_square_(5)

# control flow

values = zeros(10)

for n in 1:10
    values[n] = sin(n^2)
end

values


##
new_values = Array{Float16}(undef, 0)
for n in 1:10
    append!(new_values, sin(n^2))
end
new_values

##
f(x, y) = if x < y
    println("$x is less than $y")
elseif x > y
    println("$x is greater than $y")
else
    println("$x is equal to $y")
end

f(2, 3)

f(1, 1)

## The following sample does not work in Julia VSCode
## Julia says nn is not defined
## Maybe Jula VSCode's problem
# odds = Array{Int64}(undef, 0)
# nn = 1
# while nn <= 10
#    if isodd(nn)
#        append!(odds, nn)
#    end
#    nn = nn + 1
# end

# Random numbers

rand(5)

randn(5)

##

hist(randn(10^5))
gcf()

##
hist(randn(10^5), density=true)
gcf()

# Excercise
factorial(10)

function factorial2(x::Int)
    i = 1
    for j in 1:x
        i = i * j
    end
    return i
end

factorial2(10)

function use_while()
    odds = Array{Int64}(undef, 0)
    n = 1
    while n <= 10
        if isodd(n)
            append!(odds, n)
        end
        n = n + 1
    end
    return odds
end

use_while()

function montecarlo_pi(size_)
    pi = 0
    n = 1
    while n <= size_ 
        point = rand(2)
        if sum(point' * point) < 1
            pi = pi + 1
        end
        n = n + 1
    end
    return pi / size_ * 4
end

# 3.11
montecarlo_pi(10000)

# 3.142
montecarlo_pi(100000)