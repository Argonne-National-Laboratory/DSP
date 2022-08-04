include("../DRDCAP.jl")

generate(m,n,T,s,k) = DRDCAP(m, n, T, s, k, 0.0, "drdcap_$(m)$(n)$(T)_$(s)_$(k)")

for k in [20,50,100,200,300]
	generate(2,3,3,10,k)
	generate(2,4,3,10,k)
	generate(3,3,2,10,k)
	generate(3,4,2,10,k)
end

