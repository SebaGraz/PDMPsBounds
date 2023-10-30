
using Random 
Random.seed!(3)
d = 10
# INFECTION TIMES
tᵢ = rand(d) # epidemic values

iₒ = [1,2,4] # index obstales
tₒ = tᵢ[iₒ] .+ 0.0001 # obstacles
typeₒ = [0, 1, 0]


t = [tᵢ ; tₒ] # al times [both infection and obstacles]
f = sortperm(t) # 
b = sortperm(f) 



j = 4
println(" infection time with index ", j, " has value ", t[j])
b[j] == 1 && println("smallest time, no left neighbour");
b[j] == length(t) && println("largest time, no right neighbour");
i1, i2 = f[b[j] - 1], f[b[j] + 1]

if i1 > d
    println("The left neighbour is an obstacle of agent $(f[i1 - d]) and has value $(t[i1])")
else
    println("The left neighbour is the infection time of agent $(i1) and has value $(t[i1])")
end

if i2 > d
    println("The right neighbour is an obstacle of agent $(f[i2 - d]) and has value $(t[i2])")
   if typeₒ[i2-d] == 0
        println("and is a notification")
   else
        println("and is a removal time")
   end
else
    println("The right neighbour is the infection time of agent $(i2) and has value $(t[i2])")
end

println(" ...swapping the center with the left one")
bi1 = b[i1]
t[j], t[i1] = t[i1], t[j]
f[b[j]], f[bi1] = f[bi1], f[b[j]]
b[j], b[i1] = bi1, b[j]

println(" infection time with index ", j, " has value ", t[j])
# sanity check
fc = sortperm(t)
bc = sortperm(fc)
f == fc
bc  == b
