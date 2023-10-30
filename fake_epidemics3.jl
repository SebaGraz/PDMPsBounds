# SIMULATE EPIDEMIC AND SAVE DATA 

using Random
Random.seed!(1)
d = 20
Cij = zeros(d,d)
vicinity(i,j) = abs(i-j) < 3 ? 1/sqrt(abs(i-j)) : 0.0
for i in 1:d
    for j in 1:d
        if i == j
            Cij[i,j] = 0.0
        else
            Cij[i,j] = vicinity(i,j)*rand()*rand()
        end
    end
end
function simulate_epidemic(β, γ, Cij, I, T)
    d = length(I)
    N = zeros(d)
    R = zeros(d)
    t = 0.0
    Ξ = Vector{Tuple{Float64, Int64, Int64}}()
    for j in eachindex(I)[I .== 1]
        push!(Ξ, (t, j, 1))
    end
    while true 
        first_t = Inf
        type = 0
        j = 0
        for i in 1:d
            if I[i] == 1 
                if N[i] == 1
                    if R[i] == 1
                        continue
                    else # N -> R
                        t0 = -log(rand())/1.0
                        if t0 < first_t
                            j = i
                            type = 3
                            first_t = t0
                        end
                    end
                else # I -> N
                    t0 = -log(rand())/γ
                    if t0 < first_t
                        j = i
                        type = 2
                        first_t = t0
                    end
                end
            else # S -> I
                βi = sum(Cij[:,i])*β
                t0 = -log(rand())/βi
                if t0 < first_t
                    j = i
                    type = 1
                    first_t = t0
                end
            end        
        end
        t += first_t
        if t > T
            break
        end 
        if I[j] == 0
            I[j] = 1
        elseif N[j] == 0
            N[j] = 1
        elseif R[j] == 0
            R[j] = 1
        else
            error("")
        end
        push!(Ξ, (t, j, type))
    end
    return Ξ
end

β = 0.2
γ = 0.7
I = zeros(d)
I[5] = 1
res = simulate_epidemic(β, γ, Cij, I, T)

t1  = [res[i][1] for i in eachindex(res)[getindex.(res, 3) .== 1]]
ii = [res[i][2] for i in eachindex(res)[getindex.(res, 3) .== 1]]
ord = sortperm(ii)
t_infected = t1[ord]

t2  = [res[i][1] for i in eachindex(res)[getindex.(res, 3) .!= 1]]
ii = [res[i][2] for i in eachindex(res)[getindex.(res, 3) .!= 1]]
typeii = [res[i][3] for i in eachindex(res)[getindex.(res, 3) .!= 1]]


DATAid = Vector{Vector{Int64}}()
for i in 1:d
    obsi = findall(j -> j == i, ii)
    push!(DATAid, obsi)
end
DATA = t2

error("")


# # map index i with nothing
#             # with notification
#             # with removal 
# id_nr = Vector{Vector{Int64}}[]
# using Random, Test 
# Random.seed!(3)
# T = 4.0 # time horizon
# d = 10
# # INFECTION TIMES
# tᵢ = rand(d)*T # epidemic values
# iₒ = [1, 2, 2, 4, 0, 0] # index obstales
# tₒ = [tᵢ[1] + rand()*0.0001, tᵢ[2] + rand()*0.01,  
# tₒ = [tᵢ[iₒ[1:end-2]] .+ rand(length(iₒ)-2)*0.0001; 0.0; T + 0.0001]# obstacles
# iₒ 
# typeₒ = [0, 1, 0, 2, 2]


# # example id 3 has notification and removal times at index 1,2 and 
# # id 3 has only notification time
# # id_nr = [[],[],[1, 2],[3]]

