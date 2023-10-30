using Test
include("./gse_data.jl")
N = 199
α = 1
β = 0.001
γ = 0.15
t1, t2, id1, id2 = generate_data_gse(N, β, γ, α)
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

function next_ht(x,v,f,b)
    d = length(v)
    nt = Inf
    j = 0
    tag = 1
    for i in 1:d
        htime = ht(i, x, v, f, b)
        rtime = -log(rand()) #TODO
        if min(rtime, htime) < nt 
            if htime < rtime
                time = htime
                tag = 1
            else
                time = rtime
                tag = 2
            end
            nt, j, tag = time, i, tag
        end
    end
    return nt, j, tag
end

function next_param(p, vp)
    return -log(rand()), rand([1,2]), rand([true, false]) 
end


function swap!(i,j,x,f,b)
    x[i] = x[j]
    f[b[i]], f[b[j]] = f[b[j]], f[b[i]] 
    b[i], b[j] = b[j], b[i]
    x,f,b
end

p_swap_infections(j, k, t, f, b) = 0.5


function run_zz(p_, t_, vp_, v_, f_, b_, T, id1, id2; verbose = false)
    t, v, f, b = copy(t_), copy(v_), copy(f_), copy(b_)
    vp, p =  copy(vp_), copy(p_)
    t0 = 0.0
    d = length(v)
    println("dimension of latent space of infection times: $(d)")
    println("dimension of the parameter space $(length(p))")
    Ξ = [(t0, t_, p_),]
    while true
        δht, j, tag = next_ht(t,v,f,b)
        δrt, jp, tagp = next_param(p, vp)
        δt = min(δht, δrt)
        t[1:d] = t[1:d] + v*δt
        p = p + vp*δt
        t0 += δt
        if t0 > T
            δt = t0 - T
            t[1:d] = t[1:d] - v*δt
            t0 -= δt
            push!(Ξ, (t0, copy(t), copy(p)))
            break
        end
        push!(Ξ, (t0, copy(t), copy(p)))
        if δht < δrt # latent space
            if tag == 1 #  
                n = v[j] > 0 ? +1 : -1
                k = f[b[j] + n]
                !verbose || println("j = $(j) and k = $(k), d = $(d)")
                if k <= d
                    t[j] == t[k] || error("")
                    if rand() < p_swap_infections(j, k, t, f, b) # TODO
                        !verbose || println("swap infection times")
                        t,f,b = swap!(j, k, t,f,b)
                    else
                        v[j] *= -1
                        v[k] *= -1
                    end
                elseif k - d == d+1
                    !verbose || println("hitting 0") # check
                    v[j] == -1 || error("") # check
                    t[j] == t[k] || error("") # check
                    v[j] *= -1 
                elseif k - d == d+2
                    !verbose || println("hitting T") # check
                    v[j] == 1 || error("") # check
                    t[j] == t[k] || error("") # check
                    v[j] *= -1 
                elseif id2[k-d] == id1[j]
                    !verbose || println("Hitting its removal time, bounce off") # check
                    v[j] == 1 || error("") # check
                    v[j] *= -1
                    t[j] == t[k] || error("") # check
                elseif rand() < 0.5 # DIOSCOTINUITY, TODO
                    !verbose || println("crossing obstacle")
                    t,f,b = swap!(j, k, t,f,b)
                else
                    t[j] == t[k] || error("") # check
                    v[j] *= -1
                end
            else # tag = 2 random time latent space 
                if rand() < 0.5 #TODO
                    !verbose || println("random time latent space") # check
                    v[j] *= -1
                end
            end
        else # δrt < δht
            if !tagp 
                !verbose || println("shadow event") 
                continue
            elseif rand() < 0.5 #TODO
                !verbose || println("random time parameter space")
                vp[jp] *= -1
            end
        end
    end
    Ξ, t, v, f, b
end

T = 100.0
d = length(id1)
v = rand([+1, -1], d)
p0 = [1.0, 1.0]
vp0 = rand([-1.0, +1.0], 2)
Ξ, t, v, f, b = run_zz(p0, t0, vp0, v, f0, b0, T, id1, id2; verbose = false);
[@test t[f][i] - t[f][i-1] .>= 0 for i in 2:length(t)]

using Plots
tt = getindex.(Ξ,1)
i = 3
xxi = getindex.(getindex.(Ξ,2), i)
plot(tt,xxi)
ppi = getindex.(getindex.(Ξ,3), 1)
plot(tt,ppi)


id1[i]
i2 = findfirst(j -> j == id1[i], id2)
t2[i2]