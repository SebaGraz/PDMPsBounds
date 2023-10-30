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

function next_ht(x, v, f, b, p, vp, mp, N)
    d = length(v)
    nt = Inf
    j = 0
    ftag = 1
    for i in 1:d
        htime = ht(i, x, v, f, b)
        rtime = rand_time_tau(i, x, f, b, v, p, vp, mp, N)
        if min(rtime, htime) < nt 
            if htime < rtime
                time = htime
                tag = 1
            else
                time = rtime
                tag = 2
            end
            nt, j, ftag = time, i, tag
        end
    end
    return nt, j, ftag
end

function next_param(p, vp, t, f, b, v, mp, N, prparam)
    t1, pbs1 = rand_time_beta(p, vp, t, f, b, v, mp, N, prparam)
    t2, pbs2 = rand_time_gamma(p, vp, t, f, b, v, mp, N, prparam)
    if t1 < t2
        return t1, 1, pbs1
    else
        return t2, 2, pbs2
    end
end

function swap!(i,j,x,f,b)
    x[i] = x[j]
    f[b[i]], f[b[j]] = f[b[j]], f[b[i]] 
    b[i], b[j] = b[j], b[i]
    x,f,b
end

# Poisson time from λ(t) = max(0, -vγ*(I/(γ + vγ * t))) + max(0, vγ*sum(τcirc - τ) - vγ*sum(vτ)) 
function rand_time_gamma(p, vp, t, f, b, v, mp, N, prparam)
    λ1, ν1, λ2, ν2 =  prparam 
    I = length(v)
    γ = p[2]
    vγ = vp[2]
    τ1 = poisson_time_inv(γ, vγ, I + ν2 - 1, rand())
    aa = λ2*vγ
    for i in eachindex(v)
        aa += vγ*(t[mp[i]] - t[i])
    end
    bb = -vγ*sum(v)
    τ2 = poisson_time(aa, bb, rand())
    return min(τ1, τ2), (aa, bb, γ, vγ, I)
end
# rand_time_gamma(0.5, 1.0, ones(10), rand([-1.0,1.0]), ones(10)*2.0)


function rand_time_beta(p, vp, t, f, b, v, mp, N, prparam)
    λ1, ν1, λ2, ν2 =  prparam 
    I = length(v)
    β = p[1]
    vβ = vp[1]
    τ1 = poisson_time_inv(β, vβ, I + ν1 - 2, rand())
    aa, bb = λ1, 0.0
    for i in 1:I
        t[mp[i]] - t[i] >= 0.0 || error("tcirci - ti = $(t[mp[i]] - t[i])")
        aa += (N-I)*(t[mp[i]] - t[i])
        bb -= (N-I)*v[i]
        for j in 1:i
            if j == i
                continue
            end
            if b[i] < b[j] # τ_i < τ_j 
                if b[j] < b[mp[i]] # τ_j < τ_i^∘
                    @test t[j] - t[i] >= 0.0
                    aa += t[j] - t[i] 
                    bb += v[j] - v[i] 
                else # τ_i^circ < τ_j
                    aa += t[mp[i]] - t[i] 
                    bb -= v[i]
                end
            else # b[j] < b[i]
                if b[i] < b[mp[j]] # τ_j < τ_i^∘
                    @test t[i] - t[j] >= 0.0
                    aa += t[i] - t[j] 
                    bb += v[i] - v[j]
                else # τ_i^circ < τ_j
                    aa += t[mp[j]] - t[j]
                    bb -= v[j]
                end
            end
        end
    end
    aa *= vβ
    bb *= vβ
    τ2 = poisson_time(aa, bb, rand())
    return min(τ1, τ2), (aa, bb, β, vβ, I-1) 
end

function λp_up(δt, a, b, p, vp, c)
    return max(0, -vp*(c)/(p + vp*δt)) + max(0, a + b*δt)
end

function λp(δt, a, b, p, vp, c)
    max(0, -vp*(c/(p + vp*δt)) + a + b*δt)
end

function rand_time_tau(k, t, f, b, v, p, vp, mp, N)
    c = 0.0
    β = p[1]
    γ = p[2]
    vβ = vp[1]
    vγ = vp[2]
    for j in eachindex(v)
        if b[j] < b[k] 
            if b[k] < b[mp[j]] 
                c += 1.0
            else 
                continue
            end
        else # b[k] < b[j]
            c -= 1.0
        end
    end
    c -= N-length(v)
    aa = v[k]*(c*β - γ)
    bb = v[k]*(vβ*c - vγ) 
    τ = poisson_time(aa, bb, rand())
    return τ
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


function run_zz(p_, t_, vp_, v_, f_, b_, T, mp, N::Int64, prparam; verbose = false)
    tote = 0 
    re = 0
    t, v, f, b = copy(t_), copy(v_), copy(f_), copy(b_)
    @test sort(t) == t[f]
    vp, p =  copy(vp_), copy(p_)
    t0 = 0.0
    d = length(v)
    println("dimension of latent space of infection times: $(d)")
    println("dimension of the parameter space $(length(p))")
    Ξ = [(t0, t_, p_),]
    perc = T /10
    while true
        δht, j, tag = next_ht(t, v, f, b, p, vp, mp, N)
        δrt, jp, pbs = next_param(p, vp, t, f, b, v, mp, N, prparam)
        δrt > 0.0 || error("") 
        δht > 0.0 || error("")
        δt = min(δht, δrt)
        t[1:d] = t[1:d] + v*δt
        p = p + vp*δt
        t0 += δt
        if t0 > perc
            println("...$(round(t0/T, digits =1))...")
            perc += T/10
        end
        if t0 > T
            δt = t0 - T
            t[1:d] = t[1:d] - v*δt
            p = p - vp*δt
            t0 -= δt
            push!(Ξ, (t0, copy(t), copy(p)))
            break
        end
        push!(Ξ, (t0, copy(t), copy(p)))
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
                    error("")
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
                re += 1
                !verbose || println("random time latent space") # check
                v[j] *= -1
            end
        else # δrt < δht
            tote += 1
            prob = λp(δt, pbs...)/λp_up(δt, pbs...)
            0 <= prob <= 1.0 || error("prob = $(prob), pbs = $(pbs)") # check
            if rand() < prob  
                re += 1
                !verbose || println("random time parameter space") #check
                vp[jp] *= -1
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
    β0, γ0 = 0.01, 0.015
    # β0, γ0 = β, γ # true
    p0 = [β0, γ0]
    t1 = [t2[mpi - length(t1)] + log(rand())/γ0  for mpi in mp]
    t0 = [t1; t2] 
    f0 =  sortperm(t0)
    b0 =  sortperm(f0)
    vp1 = rand([-0.0001, 0.0001])
    vp2 = rand([-0.0001, 0.0001])
    vp0 = [vp1, vp2]
    println("check if the initial configuration of parameters and latent space is valid")
    if γ0 <= 0.0 || β0 <= 0
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
    prparam = [0.001, 1.0, 0.001, 1.0]
    Ξ, t, v, f, b = run_zz(p0, t0, vp0, v, f0, b0, T, mp, N, prparam; verbose = false);
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
pp2 = getindex.(getindex.(Ξ,3), 2)[1:100:end]
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



f2 = plot(pp1, pp2, xlabel = "β₀", ylabel = "γ₀", label = "β₀(t)-γ₀(t)", alpha = 0.5)
scatter!(f2, [βtrue],[γtrue], markershape = :circle, label = "true", color = :red)
scatter!(f2, [pp1[1]],[pp2[1]], markershape = :cross, label = "start", color = :red)

f31 = plot(tt, pp2, xlabel = "t", label = "γ₀")
hline!(f31, [γtrue], label = "γ₀ true")

f32 = plot(tt, pp1, xlabel = "t", label = "β₀")
hline!(f32, [βtrue], label = "β₀ true")

f4 = plot(tt, pp1./pp2, xlabel = "t", label = "R₀")
hline!(f4, [βtrue/γtrue], label = "R₀ true")


lf = @layout [a; b c; d e]
# lf = @layout [a b]
p = plot(f1, f31, f32, f2, f4, layout = lf, margin = 12Plots.mm, plot_title= "GSE Gibbs sampler", size = (2400, 1200))
savefig(p, "./GSE/exp2.png")
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

[t1 ; t2]