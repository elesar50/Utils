# https://www.youtube.com/watch?v=6q52CoQAt50 #To be done!
# https://github.com/dextorious/NumericalIntegration.jl

# https://github.com/JuliaDiff/ForwardDiff.jl
# For discrete points the best is:
# Integration --> NumericalIntegration.jl
# Derivation
# For continuous function then:
# Integration: QuadGK, NumericalIntegration
# Derivation: NumericalInterpolation then gradient.

using NumericalIntegration, QuadGK,  GLMakie, FiniteDiff # ReverseDiff,ForwardDiff,-->to us below to calculate the derivative on a an b points


# INTEGRATE
f(x) = sqrt(x)

a = 0.0
b = 4.0

# xs = a< b ? (a:0.01:b) : (b:0.01:a)

ya =round(f(a), digits = 2)
yb =round(f(b), digits = 2)

# dyda = round(ForwardDiff.derivative(f, a), digits = 2)
# dydb = round(ForwardDiff.derivative(f, b), digits = 2)

#To compute the area between a and b of the function f_(x)
#QuadGK integrate over a function.
area, area_error = quadgk(f, a, b)
#To integrate over points:
#
function integra(f, a, b, n)
    h = (b - a) / n

    s = (f(a) + f(b)) / 2.0
    for i in 1:n-1
        x = a + i * h
        s += f(x)
    end

    s *= h
    return s
end

n = 64

num_integral = integra(f, a, b, n)

g(x) = cos.(x)

xi = 0:0.01:4
yi = g(xi)

f1      = Figure()
ax      = Axis(f1[1,1])

lines!(xi,yi, label="Cos(x)")

Integgx = quadgk(g,a,b)
Integgx2 = integra(g,a,b,n)
#https://github.com/dextorious/NumericalIntegration.jl
# Integration by trapezoidal rule by default:
Integgy1 = integrate(xi,yi)
Integgy2 = integrate(xi,yi, SimpsonEven())
# Compute cumulative integral:
Y        = cumul_integrate(xi,yi)
lines!(xi,Y, label="Sin(x)")
YY       = cumul_integrate(xi,Y)
lines!(xi,YY, label="-Cos(x)")


#DERIVE:
# Note that what ForwardDiff package offers is not numerical differentiation. There is now a trio of approaches to differentiation:

#     symbolic: Sympy, …
#     numerical: FiniteDifferences, FiniteDiff, … ---> All need function
#     algorithmic 31: ForwardDiff, Yota, Zygote, ReverseDiff,…
#https://johnmyleswhite.github.io/FiniteDiff.jl/core/derivative/#description
#Now we derive numerically using FiniteDiff package, starting from xi and YY:
#Compute back the derivative:
# g1  = xi -> FiniteDiff.finite_difference_gradient(YY, xi)ForwardDiff,
# g1  = FiniteDiff.finite_difference_gradient(YY, xi)
# g1  = FiniteDiff.finite_difference_derivative!(YY, xi)
# g11 = g1(xi)
#Non of above permits numerical derivation...we need to interpolate the data and then applied whatever we want to get the derivative -->

# https://discourse.julialang.org/t/differentiation-without-explicit-function-np-gradient/57784

# We had YY and xi as pair of initial vector:
using Interpolations
# itp = interpolate((xi,), YY, Gridded(Linear()))
itpl= linear_interpolation(xi, YY)
# (itp(value) for each value)
# Y1  = only.(Interpolations.gradient.(Ref(itpl),3.0))
function getderivative(itp,x)
    Y = []
    for i in x
        append!(Y, only.(Interpolations.gradient.(Ref(itp),i)))
    end
    return Y
end
Y11 =getderivative(itpl,xi)
# # Y1 = Interpolations.gradient(itpl, xi)
lines!(xi,Y11, label="Sin(x)_FiniteDiff")
#And again to back to original signal:
# g2 = xi -> FiniteDiff.finite_difference_gradient(g(xi),xi)
itpll=linear_interpolation(xi,Y11)
Y111= getderivative(itpll,xi)
lines!(xi,Y111, label="Cos(x)_FiniteDiff")
f1[1,2] =Legend(f1,ax)

# THere is some numerical errors on the interpolation across the curves above on deriving it.
# Might explore other interpolation techniques but before that check forwarddiff based on the interpolated curves. To be done
# using ForwardDiff
# Y11 = ForwardDiff.derivative(itpl,3.0)
