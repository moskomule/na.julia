##

function bisection(f::Function, a, b, eps, N)
    n = 1
    p = 0.0

    while n <= N
        p = a + (b - a) / 2 # more numerically satable computation of mean

        if f(p) == 0 || abs(a - b) < eps
            return println("p is $p and the iteration number is $n")
        end

        if f(a)f(p) < 0
            b = p
        else
            a = p
        end

        n += 1
    end
    y = f(p)
    println("Did not converge. The last iteration gives $p with function value $y")
    
end

##

bisection(x -> x^5 + 2x^3 - 5x - 2, 0, 2, 10^(-4.), 20)


# obviously, if N is small, does not converge
bisection(x -> x^5 + 2x^3 - 5x - 2, 0, 2, 10^(-4.), 8)

## exercise

function bisection2(f::Function, a, b, eps, N)
    n = 1
    p = 0.0

    while n <= N
        p = a + (b - a) / 2

        if f(p) == 0 # modify stop criterion
            return println("Root found and the iteration number is $n")
        end

        if f(a)f(p) < 0
            b = p
        else
            a = p
        end

        y = f(p)
        println("At iteration $n, p=$p and f(p)=$y")

        n += 1

    end
end

bisection2(x -> x^5 + 2x^3 - 5x - 2, 0, 2, 0, 20)

## 

# newton's method

function newton(f::Function, fprime::Function, pin, eps, N)
    n = 1
    p = 0.0
    while n <= N
        p = pin - f(pin) / fprime(pin)
        if f(p) == 0 || abs(p - pin) < eps
            return println("p is $p and the iteration number is $n")
        end

        pin = p
        n += 1
    end
    y = f(p)
    println("Did not converge. The last iteration gives $p with function value of $y")
end

newton(x -> x^5 + 2x^3 - 5x - 2, x -> 5x^4 + 6x^2 - 5, 1, 10^(-4), 10)
newton(x -> x^5 + 2x^3 - 5x - 2, x -> 5x^4 + 6x^2 - 5, -2, 10^(-4), 10)

##
using PyPlot
x = range(-2, 2, length=1000)
y = map(x -> x^5 + 2x^3 - 5x - 2, x)

plot(x, y)
gcf()

##
# black-scholes formula

using Distributions

S = 7.01
K = 7.5
r = 0.0225
T = 6 / 252

stdnormal = Normal(0, 1)
phi(x) = cdf(stdnormal, x)

function c(x)
    d1 = (log(S / K) + (r + x^2 / 2) * T) / (x * sqrt(T))
    d2 = d1 - x * sqrt(T)
    return S * phi(d1) - K * exp(-r * T) * phi(d2)

end

function cprime(x)
    
    d1 = (log(S / K) + (r + x^2 / 2) * T) / (x * sqrt(T))
    d2 = d1 - x * sqrt(T)
    A = (log(S / K) + (r + x^2 / 2) * T) / (sqrt(T) * x^2)
    return S * (exp(-d1^2 / 2) / sqrt(2pi)) * (sqrt(T) - A) + K * exp(-(r * T + d2^2 / 2)) * A / sqrt(2pi)
end

newton(x -> c(x) - 0.1, x -> cprime(x), 1, 10^(-4), 50)

##
x = range(-2, 2, length=10)
x = map(x -> x^2, x)
x .+ x
x .* x

##
x = [1 2; 3 4]
x * x
x .* x # element-wise

##

# secant method

function secant(f::Function, pzero, pone, eps, N)
    n = 1
    p = 0.0

    while n <= N
        p = pone - f(pone) * (pone - pzero) / (f(pone) - f(pzero))
        if f(p) == 0 || abs(p - pone) < eps
            return println("p is $p and the iteration number is $n")
        end
        pzero = pone
        pone = p
        n += 1
    end

    y = f(p)
    println("Did not converge. The last iteration gives $p with function value $y")
    
end

secant(x -> cos(x) - x, 0.5, 1, 10^(-4), 20)

## exercise

newton(x -> log(x), x -> 1 / x, 2, 10^(-4), 20)

# x goes to negative
newton(x -> log(x), x -> 1 / x, 3, 10^(-4), 20)

##
# muller's method

function muller(f::Function, p_zero, p_one, p_two, eps, N)

    n = 1
    p = 0
    while n <= N

        c = f(p_two)
        b1 = (p_zero - p_two) * (f(p_one) - f(p_two)) / ((p_one - p_two) * (p_zero - p_one))
        b2 = (p_one - p_two) * (f(p_zero) - f(p_two)) / ((p_zero - p_two) * (p_zero - p_one))
        b = b1 - b2

        a1 = (f(p_zero) - f(p_two)) / ((p_zero - p_two) * (p_zero - p_one))
        a2 = (f(p_one) - f(p_two)) / ((p_one - p_two) * (p_zero - p_one))
        a = a1 - a2
        d = (Complex(b^2 - 4a * c))^0.5
        if abs(b - d) < abs(b + d)
            inc = 2c / (b + d)
        else
            inc = 2c / (b - d)
        end

        p = p_two - inc
        if f(p) == 0 || abs(p - p_two) < eps
            return println("p is $p and the iteration number is $n")
        end

        p_zero = p_one
        p_one = p_two
        p_two = p
        n += 1
        
    end

    y = f(p)
    println("Did not converge. The last iteration gives $p with function value $y")
end

muller(x -> x^5 + 2x^3 - 5x - 2, 0.5, 1, 1.5, 10^(-5), 10)

##

function fixed_point_iteration(g::Function, p_zero, eps, N)

    n = 1
    while n <= N
        p_one = g(p_zero)
        if abs(p_one - p_zero) < eps
            return println("p is $p_one and iteration number is $n")
        end
        p_zero = p_one
        n += 1
    end
    println("Did not converge. The last estimate is p=$p_zero")
end

fixed_point_iteration(x -> (2x^2 + 1)^(1 / 3), 1, 10^-4, 30)