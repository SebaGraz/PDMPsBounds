
using Random 
Random.seed!(3)
T = 4.0 # time horizon
d = 10
# INFECTION TIMES
tᵢ = rand(d)*T # epidemic values

iₒ = [1, 2, 4] # index obstales
tₒ = [tᵢ[iₒ] .+ rand(length(iₒ))*0.0001; 0.0; T + 0.0001]# obstacles
iₒ = [iₒ; 0; 0]
typeₒ = [0, 1, 0, 2, 2]
# 0 for notification, 1 removal, 0 boundaries

t = [tᵢ ; tₒ] # al times [both infection and obstacles]
f = sortperm(t) # 
b = sortperm(f) 

function ht(i, x, v, f, b)
    if v[i] > 0
        bj = b[i] + 1
        j = f[bj]
        if j > d
            return (x[j] - x[i])/v[i]
        elseif v[i] - v[j] > 0
            return (x[j] - x[i])/(v[i] - v[j])
        else 
            return Inf
        end
    else # v[j] < 1
        bj = b[i] - 1
        j = f[bj]
        if j > d
            return (x[j] - x[i])/v[i]
        else
            return Inf
        end
    end 
end 

function next_ht!(x,v,f,b)
    nt = Inf
    j = 0
    for i in 1:d
        htime = ht(i, x, v, f, b)
        if htime < nt 
            nt, j = htime, i
        end
    end
    return nt, j
end
function swap!(i,j,x,f,b)
    x[i] = x[j]
    f[b[i]], f[b[j]] = f[b[j]], f[b[i]] 
    b[i], b[j] = b[j], b[i]
    x,f,b
end

δt, j = next_ht!(t,v,f,b)
t[1:10] = t[1:10] + v*δt


n = v[j] > 0 ? +1 : -1
k = f[b[j] + n]
t,f,b = swap!(j, k, t,f,b)
   
[t f b]
t[f]

