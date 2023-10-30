using Distributions


function logpi(c1, c2, c3, νγ, λγ, νβ, λβ, m)
    sum(log.(c1)) - (m + νβ - 1)*log(λβ + c2) - (m + νγ)*log(λγ + c3)
end

function gse_mcmc(t0, mp, γ0, β0, νγ, λγ, νβ, λβ, m, N, iter, skip)
    Ξ = [(β0, γ0, t0), ]
    β, γ = β0, γ0
    t = copy(t0)
    ar = 0
    c1 = zeros(m-1)
    c2 = 0.0 
    c3 = 0.0
    ℓ = 1
    @show m
    _, imin = findmin(t[1:m])
    for i in 1:m
        c3 += t[mp[i]] - t[i]
        if i != imin 
            for j in 1:m
                if t[j] < t[i] < t[mp[j]]
                    c1[ℓ] += 1
                end
            end
            ℓ += 1 
        end
        for j in 1:i
            if i == j 
                continue
            else
                c2 += min(t[mp[i]], t[j]) + min(t[mp[j]], t[i]) - 2*min(t[i], t[j])  
            end
        end
        c2 += (N - m)*(t[mp[i]] - t[i])
    end
    ℓ == m  || error("ℓ equal to $ℓ")
    logπ = logpi(c1, c2, c3, νγ, λγ, νβ, λβ, m)
    for u in 1:iter-1
        β, γ, t, c1, c2, c3, logπ, ar = gse_mcmc_inner(t, mp, c1, c2, c3, logπ, m, N, νγ, λγ, νβ, λβ, ar)
        if u % skip == 0 
            push!(Ξ, (β, γ, copy(t)))
        end
    end
    println("ar is equal to $(ar/(iter-1))")
    Ξ
end

function gse_mcmc_inner(told, mp, c1old, c2old, c3old, logπold, m, N, νγ, λγ, νβ, λβ, ar)
    # update Beta
    β = rand(Gamma(m + νβ - 1, 1/(c2old + λβ)))
    # update Gamma
    γ = rand(Gamma(m + νγ, 1/(c3old + λγ)))
    # update τ_i
    k = rand(1:m)
    t = copy(told)
    ntime = -log(rand())/γ
    t[k] = t[mp[k]] - ntime
    qlogratio = γ*(told[k] - t[k])
    c1 = zero(c1old)
    c2 = 0.0 
    c3 = 0.0
    ℓ = 1
    _, imin = findmin(t[1:m])
    for i in 1:m
        c3 += t[mp[i]] - t[i]
        if i != imin    
            for j in 1:m
                if t[j] < t[i] < t[mp[j]]
                    c1[ℓ] += 1
                end
            end 
            ℓ += 1
        end
        for j in 1:i
            if i == j 
                continue
            else
                c2 += min(t[mp[i]], t[j]) + min(t[mp[j]], t[i]) - 2*min(t[i], t[j])  
            end
        end
        c2 += (N - m)*(t[mp[i]] - t[i])
    end
    ℓ == m  || error("ℓ equal to $ℓ")
    logπ = logpi(c1, c2, c3, νγ, λγ, νβ, λβ, m)
    logratio  =  logπ - logπold + qlogratio
    if log(rand()) < logratio
        ar += 1
       return β, γ, t, c1, c2, c3, logπ, ar 
    else
        return β, γ, told, c1old, c2old, c3old, logπold, ar 
    end
end


using Test
include("../GSE/gse_data.jl")
include("./../utilities.jl")
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
    data = β, γ, t1, t2, id1, id2, mp
    β0, γ0 = 0.01, 0.015
    # β0, γ0 = β, γ
    t1 = [t2[mpi - length(t1)] + log(rand())/γ0  for mpi in mp]
    t0 = [t1; t2] 
    m = length(t1)
    iter = 100000
    νγ, λγ, νβ, λβ = 1.0, 0.001, 1.0, 0.001
    res = gse_mcmc(t0, mp, γ0, β0, νγ, λγ, νβ, λβ, m, N, iter, 100)
    res, data
end
res, data = runall()
error("")
using Plots, Colors

t1 = data[3]
t2 = data[4]
mp = data[end]
dlat = length(mp)
τ1true = t1[1]
τ2true = t1[2]
τ1circ = t2[mp[1] - dlat]
τ2circ = t2[mp[2] - dlat]
ττ1 = getindex.(getindex.(res,3),1)
ττ2 = getindex.(getindex.(res,3),2)
tt = 1:length(res)

f1 = plot(tt, ττ1, label = "τ₁", color = :red, alpha =0.3, title = "latent space", xlabel = "t sampler", ylabel = "t epidemic")
hline!(f1, [t2[mp[1]-dlat]], color = :red, linestyle = :dot, label = "τ₁ circ")
hline!(f1, [τ1true], color = :red, linestyle = :dash, label = "τ₁ true")
plot!(f1, tt, ττ2, label = "τ₂", color = :blue,  alpha =0.2)
hline!(f1, [t2[mp[2]-dlat]], color = :blue, linestyle = :dot, label = "τ₂ circ")
hline!(f1, [τ2true], color = :blue, linestyle = :dash, label = "τ₂ true")




f2 = scatter(getindex.(res, 1), getindex.(res, 2),xlabel = "β₀", ylabel = "γ₀", label = "β₀(t)-γ₀(t)", alpha = 0.3, color = :blue)
scatter!(f2, [data[1]], [data[2]], markershape = :circle, label = "true", color = :red)
scatter!(f2, [getindex.(res,1)[1]], [getindex.(res,2)[1]], markershape = :cross, label = "start", color = :red)


lf = @layout [a; b c; d e]
# lf = @layout [a b]
p = plot(f1, f31, f32, f2, f4, layout = lf, margin = 12Plots.mm, plot_title= "GSE model", size = (2400, 1200), plot_title = "Metropolis withing Gibbs")
savefig(p, "./GSE/exp1.png")
error("")



f31 = plot(tt, getindex.(res,2), xlabel = "t", label = "γ₀")
hline!(f31, [data[2]], label = "γ₀ true")

f32 = plot(tt, getindex.(res,1), xlabel = "t", label = "β₀")
hline!(f32, [data[1]], label = "β₀ true")

f4 = plot(tt, getindex.(res,1)./getindex.(res,2), xlabel = "t", label = "R₀")
hline!(f4, [data[1]/data[2]], label = "R₀ true")



