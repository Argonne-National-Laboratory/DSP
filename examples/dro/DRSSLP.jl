#=
Source:
  Ntaimo, L. and S. Sen, "The 'million-variable' march for stochastic combinatorial optimization," Journal of Global Optimization, 2005.
Input:
  nJ: number of servers
  nI: number of clients
  nS: number of scenarios
Sets:
  I: clients
  J: servers
  Z: zones
Variables (1st Stage):
  x[j]: 1 if a server is located at site j, 0 otherwise
Variables (2nd Stage):
  y[i,j]: 1 if client i is served by a server at location j, 0 otherwise
  y0[j]: any overflows that are not served due to limitations in server capacity
Parameters (general):
  c[j]: cost of locating a server at location j
  q[i,j]: revenue from client i being served by server at location j
  q0[j]: overflow penalty
  d[i,j]: client i resource demand from server at location j
  u: server capacity
  v: an upper bound on the total number of servers that can be located
  w[z]: minimum number of servers to be located in zone z
  Jz[z]: the subset of server locations that belong to zone z
  p[s]: probability of occurrence for scenario s
Parameters (scenario):
  h[i,s]: 1 if client i is present in scenario s, 0 otherwise
=#

using SIPLIB
using StructJuMP
using LinearAlgebra
using Random

function write_dro(
    num_scenarios::Int, num_references::Int, 
    p::Vector{Float64}, dist::Array{Float64,2}, ϵ::Float64, 
    filename::String)

    fp = open("$filename.dro","w")
    println(fp, ϵ)
    println(fp, num_scenarios)
    for s = 1:num_scenarios
        if s > 1
            print(fp, ",")
        end
        print(fp, p[s])
    end
    print(fp, "\n")
    for i = 1:(num_scenarios+num_references)
        for s = 1:num_scenarios
            if s > 1
                print(fp, ",")
            end
            print(fp, dist[s,i])
        end
        print(fp,"\n")
    end
    close(fp)
end

function DRSSLP(nJ::Int, nI::Int, nS::Int, nK::Int, filename::String, seed::Int=1)

    Random.seed!(seed)

    J = 1:nJ
    I = 1:nI
    S = 1:nS
    K = 1:nK
    Z = []

    c = rand(40:80,nJ)
    q = rand(0:25,nI,nJ)
    q0 = ones(nJ)*1000
    d = q
    u = 1.5*sum(d)/nJ
    v = nJ
    w = NaN
    Jz = []
    h = rand(0:1,nI,nS)
    p = ones(nS) / nS
    Pr = 1.0 / (nS+nK)

    # generate many scenarios ahead
    # This allows us to consistently increase the number of scenarios for experiment.
    max_nK = Int(1e+6)
    # q_ξ = rand(0:25,nI,nJ,max_nK)
    # d_ξ = q_ξ
    # u_ξ = floor.([1.5*sum(d_ξ[:,:,s])/nJ for s=1:max_nK])
    h_ξ = rand(0:1,nI,max_nK)

    # compute Wasserstein distance
    dist = zeros(nS,nS+nK)
    for s = S, ss = S
        dist[s,ss] = sqrt(
            # norm(q[:,:,s] - q[:,:,ss])^2 
            # + norm(d[:,:,s] - d[:,:,ss])^2 
            # + norm(u[s] - u[ss])^2 
            + norm(h[:,s] - h[:,ss])^2)
    end
    for s = S, k = K
        dist[s,nS+k] = sqrt(
            # norm(q[:,:,s] - q_ξ[:,:,k])^2 
            # + norm(d[:,:,s] - d_ξ[:,:,k])^2
            # + norm(u_ξ[s] - u_ξ[k])^2  
            + norm(h[:,s] - h_ξ[:,k])^2)
    end

    write_dro(nS, nK, p, dist, 0.0001, filename)

    # construct JuMP.Model
    model = StructuredModel(num_scenarios=nS+nK)

    ## 1st stage
    @variable(model, x[j=J], Bin)
    @objective(model, Min, sum(c[j]*x[j] for j in J))
    @constraint(model, sum(x[j] for j in J) <= v)
    @constraint(model, [z=Z], sum(x[j] for j in Jz[z]) >= w[z])

    ## 2nd stage
    for ss = S
        sb = StructuredModel(parent=model, id = ss, prob = Pr)
        @variable(sb, y[i=I,j=J], Bin)
        @variable(sb, y0[j=J] >= 0)
        @objective(sb, Min, -sum(q[i,j]*y[i,j] for i in I for j in J) + sum(q0[j]*y0[j] for j in J))
        @constraint(sb, [j=J], sum(d[i,j]*y[i,j] for i in I) - y0[j] <= u*x[j])
        @constraint(sb, [i=I], sum(y[i,j] for j in J) == h[i,ss])
    end
    for k = K
        sb = StructuredModel(parent=model, id = nS+k, prob = Pr)
        @variable(sb, y[i=I,j=J], Bin)
        @variable(sb, y0[j=J] >= 0)
        @objective(sb, Min, -sum(q[i,j]*y[i,j] for i in I for j in J) + sum(q0[j]*y0[j] for j in J))
        @constraint(sb, [j=J], sum(d[i,j]*y[i,j] for i in I) - y0[j] <= u*x[j])
        @constraint(sb, [i=I], sum(y[i,j] for j in J) == h_ξ[i,k])
    end

    SIPLIB.write_smps(model, filename)

    return
end