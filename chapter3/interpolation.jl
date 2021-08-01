##

using PyPlot
using LaTeXStrings

##
function leg(x, n)
    n == 0 && return 1
    n == 1 && return x
    return ((2 * n - 1) / n) * x * leg(x, n - 1) - ((n - 1) / n) * leg(x, n - 2)
end

## plot of the first five Legendre polynomials
##

p1 = figure()
xaxis = -1:1 / 100:1
leg0 = map(x -> leg(x, 0), xaxis)
leg1 = map(x -> leg(x, 1), xaxis)
leg2 = map(x -> leg(x, 2), xaxis)
leg3 = map(x -> leg(x, 3), xaxis)
leg4 = map(x -> leg(x, 4), xaxis)

plot(xaxis, leg0, label=L"L_0(x)")
plot(xaxis, leg1, label=L"L_1(x)")
plot(xaxis, leg2, label=L"L_2(x)")
plot(xaxis, leg3, label=L"L_3(x)")
plot(xaxis, leg4, label=L"L_4(x)")
legend(loc="lower right")
gcf() 

##
close(p1)

## least-squares using Legendre polynomials

function gauss(f::Function)
    0.2369268851 * f(-0.9061798459) +
    0.2369268851 * f(0.9061798459) +
    0.5688888889 * f(0) + 
    0.4786286705 * f(0.5384693101) + 
    0.4786286705 * f(-0.5384693101)
end

##

using Base.MathConstants

gauss(x -> ((3 / 2) * x^2 - 1 / 2) * e^x)

gauss(x -> leg(x, 2) * e^x)

##

function polyLegCoeff(f::Function, n)
    global A = Array{Float64}(undef, n + 1)
    for j in 1:n + 1
        A[j] = gauss(x -> leg(x, j - 1) * f(x)) * (2 * (j - 1) + 1) / 2
    end
end

##

function polyLeg(x, n)
    sum = 0.0
    for j in 1:n + 1
        sum += A[j] * leg(x, (j - 1))
    end
    return sum
end

##

p2 = figure()
xaxis = -1:1 / 100:1
polyLegCoeff(x -> e^x,2)
deg2 = map(x -> polyLeg(x, 2), xaxis)
polyLegCoeff(x -> e^x,3)
deg3 = map(x -> polyLeg(x, 3), xaxis)
polyLegCoeff(x -> e^x,4)
deg4 = map(x -> polyLeg(x, 4), xaxis)
plot(xaxis, map(x -> e^x, xaxis), label=L"e^x")
plot(xaxis, deg2, label="Legendre least squares plot of degree 2")
plot(xaxis, deg3, label="Legendre least squares plot of degree 3")
plot(xaxis, deg4, label="Legendre least squares plot of degree 4")
legend(loc="upper left")
gcf()

##

close(p2)

##

function cheb(x, n)
    n == 0 && return 1
    n == 1 && return x
    return 2x * cheb(x, n - 1) - cheb(x, n - 2)
end

## 

p4 = figure()
cheb0 = map(x -> cheb(x, 0), xaxis)
cheb1 = map(x -> cheb(x, 1), xaxis)
cheb2 = map(x -> cheb(x, 2), xaxis)
cheb3 = map(x -> cheb(x, 3), xaxis)
cheb4 = map(x -> cheb(x, 4), xaxis)
plot(xaxis, cheb0, label=L"T_0(x)")
plot(xaxis, cheb1, label=L"T_1(x)")
plot(xaxis, cheb2, label=L"T_2(x)")
plot(xaxis, cheb3, label=L"T_3(x)")
plot(xaxis, cheb4, label=L"T_4(x)")
legend(loc="lower right")
gcf()

## 

close(p4)

##

function comp_simpson(f::Function, a, b, n)
    h = (b - a) / n
    nodes = Array{Float64}(undef, n + 1)
    for i in 1:n + 1
        nodes[i] = a + (i - 1)h
    end
    sum = f(a) + f(b)

    for i in 3:2:n - 1
        sum += 2 * f(nodes[i])
    end

    for i in 2:2:n
        sum += 4 * f(nodes[i])
    end
    return sum * h / 3
end

##

comp_simpson(x -> exp(cos(x)) * cos(x), 0, pi, 20)

## 

function polyChebCoeff(f::Function, n)
    global A = Array{Float64}(undef, n + 1)
    A[1] = comp_simpson(x -> f(cos(x)), 0, pi, 20) / pi
    for j in 2:n + 1
        A[j] = comp_simpson(x -> f(cos(x)) * cos((j - 1) * x), 0, pi, 20) * 2 / pi
    end
end

##

function polyCheb(x, n)
    sum = 0.0
    for j in 1:n + 1
        sum += A[j] * cheb(x, j - 1)
    end
    return sum
end

##

p5 = figure()
xaxis = -1:1 / 100:1
polyChebCoeff(x -> e^x,2)
deg2 = map(x -> polyCheb(x, 2), xaxis)
polyChebCoeff(x -> e^x,3)
deg3 = map(x -> polyCheb(x, 3), xaxis)
polyChebCoeff(x -> e^x,4)
deg4 = map(x -> polyCheb(x, 4), xaxis)
plot(xaxis, map(x -> e^x, xaxis), label=L"e^x")
plot(xaxis, deg2, label="Chebyshev least squares plot of degree 2")
plot(xaxis, deg3, label="Chebyshev least squares plot of degree 3")
plot(xaxis, deg4, label="Chebyshev least squares plot of degree 4")
legend(loc="upper left")
gcf()

##

close(p5)

