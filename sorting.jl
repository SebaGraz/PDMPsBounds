using Random
Random.seed!(1)
d = 5
x = [0.24, 0.35, 0.31, 0.01, 0.49] # round.(rand(d), digits=2)
f = sortperm(x)
b = sortperm(f)

i = 1
println(x, " ", x[i])
println(sort(x))
println(f) # locate n'th 
println(b) # locate order [inverse funciton]

f1, f2 = f[b[i]-1], f[b[i]+1]
println(i, "(",x[i], ")", " is between ", f1,"(",x[f1],") and ", f2, "(", x[f2],")")
x[f]
x[f[b[1]]]

i = 1 # index i
x[i] # with value xi
b[i] # is b[i] smallest value
x[f[b[i]]] # and is equal to  x[f[b[i]]] 

println("swapping with the right-neighbour")
bj = b[i] + 1
j = f[bj]
y = x[i]
x[i], x[j] = x[j], x[i]
f[b[i]], f[bj] = f[bj], f[b[i]]
b[i], b[j] = b[j], b[i]


