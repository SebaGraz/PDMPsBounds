
"""
poisson_time_inv(c, x, v, u)


Obtaining waiting time for inhomogeneous Poisson Process
    with rate of the form λ(t) = -v(c/(x + vt))^+, `x + vt > 0`,`c>0`, `u` uniform random variable
    noitce that the rate goes to infinity as x + vt -> 0.0. 
"""
function poisson_time_inv(x, v, c, u)
    if v > 0.0
        return Inf
    else
        return (exp(log(u)/c)*x - x)/v
    end
end



"""
poisson_time_inv2(c, a, b, u)


Obtaining waiting time for inhomogeneous Poisson Process
    with rate of the form λ(t) = (1/(a + bt))^+, 
"""
function poisson_time_inv2(a, b, u)
    e = -log(u)
    if a<0 || abs(b) < eps()
        return Inf
    else
        res = (exp(e*b + log(a)) - a)/b
        if res < 0
            error("res = $(res), a = $(a), b = $b")
        end
        return res
    end
end



# x = 1/2
# v = -1.0
# c = 2.0
# lambda_1(x,v, c) = max(-c*v/x, 0.0)
# using Plots
# xx = 0.001:0.001:1.0
# plot(xx, [lambda_1(x,v,c) for x in xx])
# t = poisson_time_inv(x, v, c, rand())
# x = x + v*t
# vline!([x])


"""
poisson_time(a, b, u)

Obtaining waiting time for inhomogeneous Poisson Process
with rate of the form λ(t) = (a + b*t)^+, `a`,`b` ∈ R, `u` uniform random variable
"""
function poisson_time(a, b, u)
if b > 0
    if a < 0
        return sqrt(-log(u)*2.0/b) - a/b
    else # a[i]>0
        return sqrt((a/b)^2 - log(u)*2.0/b) - a/b
    end
elseif b == 0
    if a > 0
        return -log(u)/a
    else # a[i] <= 0
        return Inf
    end
else # b[i] < 0
    if a <= 0
        return Inf
    elseif -log(u) <= -a^2/b + a^2/(2*b)
        return -sqrt((a/b)^2 - log(u)*2.0/b) - a/b
    else
        return Inf
    end
end
end

"""
poisson_time((a, b, c), u)

Obtaining waiting time for inhomogeneous Poisson Process
with rate of the form λ(t) = c + (a + b*t)^+, 
where `c`> 0 ,`a, b` ∈ R, `u` uniform random variable
"""
function poisson_time((a,b,c)::NTuple{3}, u)
    if b > 0
        if a < 0
            if -c*a/b + log(u) < 0.0
                return sqrt(-2*b*log(u) + c^2 + 2*a*c)/b - (a+c)/b
            else
                return -log(u)/c
            end
        else # a >0
            return sqrt(-log(u)*2.0*b + (a+c)^2)/b - (a+c)/b
        end
    elseif b == 0
        if a > 0
            return -log(u)/(a + c)
        else # a <= 0
            return -log(u)/c
        end
    else # b < 0
        if a <= 0.0 
            return -log(u)/c
        elseif  - c*a/b - a^2/(2*b)  + log(u) > 0.0
            return +sqrt((a+c)^2 - 2.0*log(u)*b)/b - (a+c)/b
        else
            return (-log(u)+ a^2/(2*b))/c
        end
    end
end