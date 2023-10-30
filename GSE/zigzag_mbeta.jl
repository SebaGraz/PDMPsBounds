using Test
include("./gse_data.jl")
include("./../utilities.jl")


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

function next_ht(x, v, f, b, γ, vγ, mp, N, prparam)
    # return Inf, 0, 1, ()
    d = length(v)
    nt = Inf
    j = 0
    ftag = 1
    fpbs = (0.0, 0.0, 0.0, 0.0)
    for i in 1:d
        htime = ht(i, x, v, f, b)
        rtime, rpbs = rand_time_tau(i, x, f, b, v, γ, vγ, mp, N, prparam)
        if min(rtime, htime) < nt 
            if htime < rtime
                time = htime
                tag = 1
            else
                time = rtime
                tag = 2                
            end
            pbs = rpbs
            nt, j, ftag, fpbs = time, i, tag, pbs
        end
    end
    return nt, j, ftag, fpbs
end

function next_param(γ::Float64, vγ::Float64, t, f, b, v, mp, N, prparam::Tuple{Float64, Float64, Float64, Float64})
    t1, pbs1 = rand_time_gamma(γ, vγ, t, f, b, v, mp, N, prparam)
    return t1, pbs1
end

function swap!(i,j,x,f,b)
    x[i] = x[j]
    f[b[i]], f[b[j]] = f[b[j]], f[b[i]] 
    b[i], b[j] = b[j], b[i]
    x,f,b
end

function rand_time_gamma(γ::Float64, vγ::Float64, t, f, b, v, mp, N, prparam::Tuple{Float64, Float64, Float64, Float64})
    λ1, ν1, λ2, ν2 =  prparam 
    I = length(v)
    τ1 = poisson_time_inv(γ, vγ, I + ν2 - 1, rand())
    aa = λ2*vγ
    for i in eachindex(v)
        aa += vγ*(t[mp[i]] - t[i])
    end
    bb = -vγ*sum(v)
    τ2 = poisson_time(aa, bb, rand())
    return min(τ1, τ2), (aa, bb, γ, vγ, I + ν2 - 1)
end
# rand_time_gamma(0.5, 1.0, ones(10), rand([-1.0,1.0]), ones(10)*2.0)


function λγ_up(δt, a, b, γ, vγ, c)
    return max(0, -vγ*(c)/(γ + vγ*δt)) + max(0, a + b*δt)
end

function λγ(δt, a, b, γ, vγ, c)
    max(0, -vγ*(c/(γ + vγ*δt)) + a + b*δt)
end



function rand_time_tau(k, t, f, b, v, γ, vγ, mp, N, prparam)
    λ1, ν1, λ2, ν2 =  prparam 
    m = length(v)
    c = 0.0
    aa = λ1
    bb = 0.0
    for i in eachindex(v)
        if b[i] < b[k] 
            if b[k] < b[mp[i]] 
                c += 1.0
            else 
                continue
            end
        else # b[k] < b[j]
            c -= 1.0
        end
        for j in 1:i
            if i != j
                if b[i] < b[j]
                    if b[j] < b[mp[i]]
                        aa += t[j] - t[i]
                        bb += v[j] - v[i]
                    else # b[mp[i]] < b[j]
                        aa += t[mp[i]] - t[i]
                        bb -= v[i]
                    end
                else # b[j] < b[i] 
                    if b[i] < b[mp[j]]
                        aa += t[i] - t[j]
                        bb += v[i] - v[j]
                    else # b[mp[j]] < b[i]
                        aa += t[mp[j]] - t[j]
                        bb -= v[j] 
                    end
                end
            end
        end
        aa += (N-m)*(t[mp[i]] - t[i])
        bb -= (N-m)*(v[i])
    end
    c -= N-m
    c *= (ν1 + m - 1)*v[k] 
    aa = aa/c
    bb = bb/c
    τ1 = poisson_time_inv2(aa, bb, rand())
    @test τ1 > 0.0
    τ2 = poisson_time(-γ*v[k], -vγ*v[k], rand())
    τ = min(τ1, τ2) 
    return τ, (aa, bb, -γ*v[k], -vγ*v[k])
end

function λtau_up(δt, aa, bb, cc, dd)
    return max(1/(aa + bb*δt), 0) + max(cc + dd*δt, 0.0)
end

function λtau(δt, aa, bb, cc, dd)
    max(1/(aa + bb*δt) + cc + dd*δt, 0.0)
end

function p_swap_infections(j, k, t, f, b, v, mp) 
    if k <= length(v) # swap infection times (homogeneous mixing happens always)
        return 1.0
    elseif mp[j] == k # collision with its own removal time
        return 1.0
    else # collision with removal time
        c = 0.0
        for i in eachindex(v)
            if b[i] < b[j] < b[mp[i]]
                c += 1 
            end 
        end
        if v[j] > 0.0
            return (c-1)/c
        else
            return (c+1)/c
        end
    end
end


function run_zz(γ::Float64, t_, vγ::Float64, v_, f_, b_, T, mp, N::Int64, prparam::Tuple{Float64, Float64, Float64, Float64}; verbose = false)
    tote = 0 
    re = 0
    t, v, f, b = copy(t_), copy(v_), copy(f_), copy(b_)
    @test sort(t) == t[f]
    t0 = 0.0
    d = length(v)
    println("dimension of latent space of infection times: $(d)")
    println("marignalizing over β")
    Ξ = [(t0, t_, γ),]
    perc = T /10
    while true
        δht, j, tag, pbs1 = next_ht(t, v, f, b, γ, vγ, mp, N, prparam)
        δrt, pbs2 = next_param(γ, vγ, t, f, b, v, mp, N, prparam)
        δrt > 0.0 || error("") 
        δht > 0.0 || error(" δht = $(δht)")
        δt = min(δht, δrt)
        t[1:d] = t[1:d] + v*δt
        γ += vγ*δt
        t0 += δt
        if t0 > perc
            println("...$(round(t0/T, digits =1))...")
            perc += T/10
        end
        if t0 > T
            δt = t0 - T
            t[1:d] = t[1:d] - v*δt
            γ = γ - vp*δt
            t0 -= δt
            push!(Ξ, (t0, copy(t), γ))
            break
        end
        push!(Ξ, (t0, copy(t), γ))
        if δht < δrt # latent space
            if tag == 1 # hitting time 
                n = v[j] > 0 ? +1 : -1
                k = f[b[j] + n]
                !verbose || println("j = $(j) and k = $(k), d = $(d)")
                if k <= d
                    t[j] = t[k] # || error("t[j] == t[k]; $(t[j]) == $(t[k])")
                    if rand() < p_swap_infections(j, k, t, f, b, v, mp)  
                        !verbose || println("swap infection times")
                        t,f,b = swap!(j, k, t,f,b)
                    else
                        v[j] *= -1
                        v[k] *= -1
                    end
                elseif k - d == d+1
                    !verbose || println("hitting lower bound $(t[k])") # check
                    v[j] == -1 || error("") # check
                    t[j] == t[k] || error("") # check
                    v[j] *= -1 
                elseif k - d == d+2
                    error("don't need the boundary")
                    !verbose || println("hitting T") # check
                    v[j] == 1 || error("") # check
                    t[j] == t[k] || error("") # check
                    v[j] *= -1 
                elseif mp[j] == k 
                    !verbose || println("Hitting its removal time, bounce off") # check
                    v[j] == 1 || error("") # check
                    v[j] *= -1
                    t[j] == t[k] || error("") # check
                elseif rand() < p_swap_infections(j, k, t, f, b, v, mp) # DIOSCOTINUITY, TODO
                    !verbose || println("crossing obstacle")
                    t,f,b = swap!(j, k, t,f,b)
                else
                    t[j] == t[k] || error("") # check
                    v[j] *= -1
                end
            else # tag = 2 random time latent space 
                tote += 1
                prob = λtau(δt, pbs1...)/λtau_up(δt, pbs1...)
                0 <= prob <= 1.00001 || error("prob $(prob)") # check
                if rand() < prob  
                    re += 1
                    !verbose || println("random time latent space") # check
                    v[j] *= -1
                else
                    !verbose || println("shadow event") 
                    continue
                end
            end
        else # δrt < δht
            tote += 1
            prob = λγ(δt, pbs2...)/λγ_up(δt, pbs2...)
            0 <= prob <= 1.0 || error("") # check
            if rand() < prob  
                re += 1
                !verbose || println("random time parameter space") #check
                vγ *= -1
            else
                !verbose || println("shadow event") 
                continue
            end
        end
    end
    println("ratio between random events and total random times: $(re/tote)")
    Ξ, t, v, f, b
end


function runall()
    N = 200
    β =  0.001
    γ =  0.15
    t1, t2, id1, id2 = generate_data_gse(N, β, γ) 
    mp = Vector{Int64}()
    for i in eachindex(t1)
        id = findfirst(j -> j == id1[i], id2)
        push!(mp, id + length(t1))
        t1[i] < t2[id] || error("")
    end
    data = β, γ, t1, t2, id1, id2, mp
    # ttrue = [t1 ; t2] # all times [both infection and obstacles]
    # ftrue = sortperm(ttrue) #
    # btrue = sortperm(ftrue) 

    T = 1000.0
    d = length(id1)
    v = rand([-1.0, 1.0], d)
    γ0 = 0.015
    β0, γ0 = β, γ # true
    # t1 = [t2[mpi - length(t1)] + log(rand())/γ0  for mpi in mp]
    t0 = [t1; t2] 
    f0 =  sortperm(t0)
    b0 =  sortperm(f0)
    vγ0 = rand([-0.01, 0.01])
    println("check if the initial configuration of parameters and latent space is valid")
    if γ0 <= 0.0 
        error("parameters out of the domain")
    end
    ni = 1
    inf_times = t1
    dim = length(t1)
    rem_times = t2[1:end-2]
    tall = [inf_times ; rem_times]
    fall = sortperm(tall)
    ni = 0
    for i in eachindex(fall)
        if fall[i] > dim
            ni -= 1
        else
            ni += 1
        end
        if ni <= 0
            if length(fall) == i
                continue
            else
                error("Likelihood at the starting value equal to 0: this is not a valid configuration...")
            end
        end 
    end
    λ1, λ2 = 0.001, 0.001
    ν1, ν2 = 1.0, 1.0
    par_prior = λ1, ν1, λ2, ν2
    Ξ, t, v, f, b = run_zz(γ0, t0, vγ0, v, f0, b0, T, mp, N,  par_prior; verbose = false);
    [@test t[f][i] - t[f][i-1] .>= 0 for i in 2:length(t)]
    return Ξ, t, v, f, b, data
end
Ξ, t, v, f, b, data = runall()

using Plots  
mp = data[end]
dlat = length(mp)
t1 = data[3]
d = length(t1)
t2 = data[4]
tt = getindex.(Ξ,1)[1:100:end]
pp1 = getindex.(getindex.(Ξ,3), 1)[1:100:end]
γtrue = data[2]
βtrue = data[1]

ττ1 = getindex.(getindex.(Ξ,2), 1)[1:100:end]
ττ2 = getindex.(getindex.(Ξ,2), 2)[1:100:end]
τ1true = t1[1]
τ2true = t1[2]
τ1circ = t2[mp[1] - dlat]
τ2circ = t2[mp[2] - dlat]

f1 = plot(tt, ττ1, label = "τ₁", color = :red, alpha =0.3, title = "latent space", xlabel = "t sampler", ylabel = "t epidemic")
hline!(f1, [t2[mp[1]-dlat]], color = :red, linestyle = :dot, label = "τ₁ circ")
hline!(f1, [τ1true], color = :red, linestyle = :dash, label = "τ₁ true")
plot!(f1, tt, ττ2, label = "τ₂", color = :blue,  alpha =0.2)
hline!(f1, [t2[mp[2]-dlat]], color = :blue, linestyle = :dot, label = "τ₂ circ")
hline!(f1, [τ2true], color = :blue, linestyle = :dash, label = "τ₂ true")


f31 = plot(tt, pp1, xlabel = "t", label = "γ₀")
hline!(f31, [γtrue], label = "γ₀ true")


lf = @layout [a; b]
# lf = @layout [a b]
p = plot(f1, f31, layout = lf, margin = 12Plots.mm, plot_title= "GSE Gibbs sampler", size = (2400, 1200))
savefig(p, "./GSE/exp_mbeta.png")
error("")


inf_times = data[3]
rem_times = data[4][1:end-2]
dim = length(inf_times)
tall = [inf_times ; rem_times]
f = sortperm(tall)
inf_t = zeros(dim*2 + 1)
succ_t = fill(200, dim*2 + 1)
ni = 0
for i in eachindex(f)
    if f[i] > dim
        ni -= 1
        inf_t[i+1] = inf_t[i] - 1
        succ_t[i+1] = succ_t[i] + 1
    else
        ni += 1
        inf_t[i+1] = inf_t[i] + 1
        succ_t[i+1] = succ_t[i] - 1
    end
    if ni <= 0
        if length(f) == i
            continue
        else
            error("iteration $i")
        end
    end 
end
ni
scatter(tall[f], inf_t[2:end], xlabel = "t", ylabel = "|I(t)|")

# scatter(tall[f], succ_t, legend = false)
