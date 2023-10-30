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

function next_param(p, vp, t, f, b, v, mp, N)
    t1, pbs1 = rand_time_beta(p, vp, t, f, b, v, mp, N)
    t2, pbs2 = rand_time_gamma(p, vp, t, f, b, v, mp, N)
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
function rand_time_gamma(p, vp, t, f, b, v, mp, N)
    I = length(v)
    τ1 = poisson_time_inv(p[2], vp[2], I, rand())
    aa = 0.0
    for i in eachindex(v)
        aa += vp[2]*(t[mp[i]] - t[i])
    end
    bb = -vp[2]*sum(v)
    τ2 = poisson_time(aa, bb, rand())
    return min(τ1, τ2), (aa, bb, p[2], vp[2], I)
end
# rand_time_gamma(0.5, 1.0, ones(10), rand([-1.0,1.0]), ones(10)*2.0)


function rand_time_beta(p, vp, t, f, b, v, mp, N)
    I = length(v)
    τ1 = poisson_time_inv(p[1], vp[1], I-1, rand())
    aa, bb = 0.0, 0.0
    for i in 1:I
        t[mp[i]] - t[i] >= 0.0 || error("tcirci - ti = $(t[mp[i]] - t[i])")
        aa += (N-I)*(t[mp[i]] - t[i])
        bb -= (N-I)*v[i]
        for j in 1:I
            if j == i
                continue
            end
            if b[i] < b[j] # τ_i < τ_j 
                if b[j] < b[mp[i]] # τ_j < τ_i^∘
                    aa += t[j] - t[i] 
                    bb += v[j] - v[i] 
                else # τ_i^circ < τ_j
                    aa += t[mp[i]] - t[i] 
                    bb -= v[i]
                end
            else # b[j] < b[i]
                if b[i] < b[mp[j]] # τ_j < τ_i^∘
                    aa += t[i] - t[j] 
                    bb += v[i] - v[j]
                else # τ_i^circ < τ_j
                    aa += t[mp[j]] - t[j]
                    bb -= v[j]
                end
            end
        end
    end
    aa *= vp[1]
    bb *= vp[1]
    τ2 = poisson_time(aa, bb, rand())
    return min(τ1, τ2), (aa, bb, p[1], vp[1], I-1) 
end

function λp_up(δt, a, b, p, vp, c)
    return max(0, -vp*(c)/(p + vp*δt)) + max(0, a + b*δt)
end

function λp(δt, a, b, p, vp, c)
    pnew = p + vp*δt
    max(0, -vp*(c/pnew) + a)
end
function rand_time_tau(k, t, f, b, v, p, vp, N)
    c = 0.0
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
    aa = v[k]*(c*p[1] - p[2])
    bb = v[k]*(vp[1]*c - vp[2]) 
    τ2 = poisson_time(aa, bb, rand())
    return min(τ1, τ2)
end


function p_swap_infections(j, k, t, f, b, v, mp) 
    if k <= length(v) # swap infection times (homogeneous mixing happens always)
        return 1.0
    elseif mp[j] == k # collision with its own removal time
        return 1.0
    else # collision with removal time
        c = 0.0
        for i in eachindex(v)
            if b[i] < b[k] < b[mp[i]]
                c += 1 
            end 
        end
        if v[j] > 0.0
            return 1.0
        else
            return 1.0
        end
    end
end


function run_zz(p_, t_, vp_, v_, f_, b_, T, mp, N::Int64; verbose = false)
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
        δht, j, tag = next_ht(t,v,f,b)
        δrt, jp, pbs = next_param(p, vp, t, f, b, v, mp, N)
        δrt >= 0.0 || error("") 
        if δht < 0.0 
            for i in eachindex(b)
                @test sort(t)[i] == t[f][i]
            end
            error("$(δht)")
        end
        δt = min(δht, δrt)
        @show δrt
        @show δht
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
                !verbose || println("random time latent space") # check
                v[j] *= -1
            end
        else # δrt < δht
            prob = λp(δt, pbs...)/λp_up(δt, pbs...)
            prob < 1.0 || error("") # check
            if rand() < prob  
                !verbose || println("random time parameter space") #check
                vp[jp] *= -1
            else
                !verbose || println("shadow event") 
                continue
            end
        end
    end
    Ξ, t, v, f, b
end


function runall()
    N = 200
    β = 0.001
    γ = 0.15
    t1, t2, id1, id2 = generate_data_gse(N, β, γ) 
    mp = Vector{Int64}()
    for i in eachindex(t1)
        id = findfirst(j -> j == id1[i], id2)
        push!(mp, id + length(t1))
        t1[i] < t2[id] || error("")
    end
    data = t1, t2, id1, id2, mp
    # ttrue = [t1 ; t2] # all times [both infection and obstacles]
    # ftrue = sortperm(ttrue) #
    # btrue = sortperm(ftrue) 

    T = 100.0
    d = length(id1)
    v = rand([-1.0, 1.0], d)
    # β0, γ0 = 0.01, 1.5
    β0, γ0 = β, γ
    p0 = [β0, γ0]
    # t1 = [t2[mpi - length(t1)] + log(rand())/γ0  for mpi in mp]
    t0 = [t1; t2]
    f0 =  sortperm(t0)
    b0 =  sortperm(f0)
    

    vp1 = rand([-0.001, 0.001])
    vp2 = rand([-0.1, 0.1])
    vp0 = [vp1, vp2]
    Ξ, t, v, f, b = run_zz(p0, t0, vp0, v, f0, b0, T, mp, N; verbose = false);
    [@test t[f][i] - t[f][i-1] .>= 0 for i in 2:length(t)]
    return Ξ, t, v, f, b, data
end
Ξ, t, v, f, b, data = runall()
error("")




using Plots
mp = data[end]
t1 = data[1]
d = length(t1)
t2 = data[2]
tt = getindex.(Ξ,1)
i = 1 
f1 = plot(tt, getindex.(getindex.(Ξ,2), 1), label = "τ₁", color = :red )
hline!(tt, [t2[mp[1]-71]], color = :red)
plot!(f1, tt, getindex.(getindex.(Ξ,2), 2), label = "τ₂", color = :blue)


pp1 = getindex.(getindex.(Ξ,3), 1)
f2 = plot(tt,pp1)
hline!([0.001])
pp2 = getindex.(getindex.(Ξ,3), 2)
f3 = plot(tt,pp2)
hline!([0.15])
inf_times = data[1]
rem_times = data[2][1:end-2]
dim = length(inf_times)
tall = [inf_times ; rem_times]
f = sortperm(tall)
inf_t = zero(tall)
inf_t[1] = 1.0
succ_t = fill(200 - 1, d)
for i in eachindex(f)[2:end]
    if i == 1
        f[1] == 1 || error("")
    elseif f[i] > dim
        inf_t[i] = inf_t[i-1] - 1
        succ_t[i] = succ_t[i-1] + 1
    else
        inf_t[i] = inf_t[i-1] + 1
        succ_t[i] = succ_t[i-1] - 1
    end
end
scatter(tall[f], inf_t, xlabel = "t", ylabel = "|I(t)|")

# scatter(tall[f], succ_t, legend = false)

