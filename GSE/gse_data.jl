# See https://link.springer.com/article/10.1007/s11222-005-4074-7
using Random


function simulate_epidemic(β, γ, N; verbose = true)
    println("population size = $(N )")
    I = zeros(N)
    I[1] = 1
    println("infected initial population = $(sum(I))")
    d = length(I)
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
                if R[i] == 1
                    continue
                else # I -> R
                    t0 = -log(rand())/γ
                    if t0 < first_t
                        j = i
                        type = 2
                        first_t = t0
                    end
                end
            else # S -> I
                βi = (sum(I) - sum(R))*β
                 t0 = -log(rand())/βi
                if t0 < first_t
                    j = i
                    type = 1
                    first_t = t0
                end
            end        
        end
        t += first_t
        if I[j] == 0
            I[j] = 1
        elseif R[j] == 0
            R[j] = 1
        else
            error("")
        end
        push!(Ξ, (t, j, type))
        if (sum(I) - sum(R)) == 0
            println("last removed individual: $(t)")
            println("total number of infected $(sum(I))")
            break
        end
    end
    return Ξ
end

function generate_data_gse(N, β, γ)
    Random.seed!(1)
    res = simulate_epidemic(β, γ, N)
    t1  = [res[i][1] for i in eachindex(res)[getindex.(res, 3) .== 1]]
    id1 = [res[i][2] for i in eachindex(res)[getindex.(res, 3) .== 1]]
    # ii = [res[i][2] for i in eachindex(res)[getindex.(res, 3) .== 1]]
    # ord = sortperm(ii)
    # t_infected = t1[ord]
    
    t2  = [res[i][1] for i in eachindex(res)[getindex.(res, 3) .!= 1]]
    push!(t2, -Inf)
    push!(t2, maximum(t2))
    id2 = [res[i][2] for i in eachindex(res)[getindex.(res, 3) .!= 1]]
    # typeii = [res[i][3] for i in eachindex(res)[getindex.(res, 3) .!= 1]]
    # DATAid = Vector{Vector{Int64}}()
    # for i in 1:length(t1)
    #     obsi = findall(j -> j == i, ii)
    #     push!(DATAid, obsi)
    # end
    return t1, t2, id1, id2
end

# N = 200
# β = 0.001
# γ = 0.15
# t1, t2, DATAid = generate_data_gse(N, β, γ)
