using Test
include("./gse_data.jl")
N = 199
α = 1
β = 0.001
γ = 0.15
t1, t2, DATAid = generate_data_gse(N, β, γ, α)
t0 = [t1 ; t2] # all times [both infection and obstacles]
f0 = sortperm(t0) # 
b0 = sortperm(f0) 

function ht(i, x, v, f, b)
    d = length(v)
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
    d = length(v)
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


function run_Its(t, v, f, b, T, id; verbose = true)
    t0 = 0.0
    d = length(v)
    while true
        δt, j = next_ht!(t,v,f,b)
        t[1:d] = t[1:d] + v*δt
        n = v[j] > 0 ? +1 : -1
        k = f[b[j] + n]
        !verbose || println("j = $(j) and k = $(k), d = $(d)")
        if k <= d
            !verbose || println("swap infection times")
            t,f,b = swap!(j, k, t,f,b)
        elseif k - d == d+1
            !verbose || println("hitting 0")
            v[j] == -1 || error("")
            v[j] *= -1 
        elseif k - d == d+2
            !verbose || println("hitting T")
            v[j] == 1 || error("")
            v[j] *= -1 
        elseif id[k-d] == j
            !verbose || println("Hitting removal time, bounce off")
            v[j] == 1 || error("")
            v[j] *= -1
        else
            !verbose || println("crossing obstacle")
            t,f,b = swap!(j, k, t,f,b)
        end
        t0 += δt
        if t0 > T
            break
        end
    end
    t, v, f, b
end

T = 10.0
d = length(DATAid)
v = rand([+1, -1], d)
t, v, f, b = run_Its(t0, v, f0, b0, T, DATAid)
[@test t[f][i] - t[f][i-1] .>= 0 for i in 2:length(t)]
