using Plots
using DelimitedFiles
using Printf

data = readdlm("tmp.csv", ',', Float64, '\n', header=true)

datadict = Dict{String,Dict}()
iter = Dict{Vector{Int},Int}()
UB = Dict{Vector{Int},Float64}()
LB = Dict{Vector{Int},Float64}()
time = Dict{Vector{Int},Float64}()
datadict["iter"] = iter
datadict["UB"] = UB
datadict["LB"] = LB
datadict["time"] = time
for i = 1:size(data[1],1)
    datadict["iter"][convert(Vector{Int},data[1][i,[1,2,3]])] = Int(data[1][i,4])
    datadict["UB"][convert(Vector{Int},data[1][i,[1,2,3]])] = data[1][i,5]
    datadict["LB"][convert(Vector{Int},data[1][i,[1,2,3]])] = data[1][i,6]
    datadict["time"][convert(Vector{Int},data[1][i,[1,2,3]])] = data[1][i,7]
end

###################
## Generate figures
###################

instances = [233,243,332,342]
K = [20,50,100,200,300]
eps = [1,100,500,1000]

for instance in instances
    vals = zeros(length(eps),length(K))
    for i = 1:length(K), j = 1:length(eps)
        vals[j,i] = datadict["LB"][[instance,K[i],eps[j]]]
    end
    fs = 14
    plot(eps, vals, 
        tickfontsize=fs, legendfontsize=fs, guidefontsize=fs,
        label=["k=$k" for k=K'], legend=:topleft,
        xlabel="Epsilon", ylabel="Lower bound")
    savefig("dcap$(instance)_LB.pdf")
end

for instance in instances
    vals = zeros(length(eps),length(K))
    for i = 1:length(K), j = 1:length(eps)
        vals[j,i] = datadict["time"][[instance,K[i],eps[j]]]
    end
    fs = 14
    plot(eps, vals, 
        tickfontsize=fs, legendfontsize=fs, guidefontsize=fs,
        label=["k=$k" for k=K'], legend=:topright,
        xlabel="Epsilon", ylabel="Solution time (s)")
    savefig("dcap$(instance)_time.pdf")
end

#######################
## Generate Latex table
#######################

println("\\begin{tabular}{lrrrrr}")
println("\\hline")
println("Name & Iteration & UB & LB & Gap (\\%) & Time (s) \\\\")
for instance in [233,243]
    for k in [20,50,100,200,300]
        println("\\hline")
        for e in [1,100,500,1000]
            iter = datadict["iter"][[instance,k,e]]
            ub = datadict["UB"][[instance,k,e]]
            lb = datadict["LB"][[instance,k,e]]
            gap = (ub-lb)/abs(1e-10+ub)*100
            time = datadict["time"][[instance,k,e]]
            @printf("dcap%d\\_%d\\_%d & %d & %.2f & %.2f & %.2f & %d \\\\\n", instance, k, e, iter, ub, lb, gap, time)
        end
    end
end
println("\\hline")
println("\\end{tabular}")