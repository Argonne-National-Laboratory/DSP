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

function DRSSLP(nJ::Int, nI::Int, nS::Int, nK::Int, filename::String, 
  relax_xbin::Bool = false, relax_ybin::Bool = false, seed::Int=1)

    Random.seed!(seed)

    J = 1:nJ
    I = 1:nI
    S = 1:nS # number of references
    K = 1:nK # number of discretization points
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
    h_ξ = rand(0:1,nI,max_nK)

    # construct JuMP.Model
    model = StructuredModel(num_scenarios=nS+nK)

    ## 1st stage
    if relax_xbin
      @variable(model, 0 <= x[j=J] <= 1)
    else
      @variable(model, x[j=J], Bin)
    end
    @objective(model, Min, sum(c[j]*x[j] for j in J))
    @constraint(model, sum(x[j] for j in J) <= v)
    @constraint(model, [z=Z], sum(x[j] for j in Jz[z]) >= w[z])

    ## 2nd stage
    for ss = S
        sb = StructuredModel(parent=model, id = ss, prob = Pr)
        if relax_ybin
          @variable(sb, 0 <= y[i=I,j=J] <= 1)
        else
          @variable(sb, y[i=I,j=J], Bin)
        end
        @variable(sb, y0[j=J] >= 0)
        @objective(sb, Min, -sum(q[i,j]*y[i,j] for i in I for j in J) + sum(q0[j]*y0[j] for j in J))
        @constraint(sb, [j=J], sum(d[i,j]*y[i,j] for i in I) - y0[j] <= u*x[j])
        @constraint(sb, [i=I], sum(y[i,j] for j in J) == h[i,ss])
    end
    for k = K
        sb = StructuredModel(parent=model, id = nS+k, prob = Pr)
        if relax_ybin
          @variable(sb, 0 <= y[i=I,j=J] <= 1)
        else
          @variable(sb, y[i=I,j=J], Bin)
        end
        @variable(sb, y0[j=J] >= 0)
        @objective(sb, Min, -sum(q[i,j]*y[i,j] for i in I for j in J) + sum(q0[j]*y0[j] for j in J))
        @constraint(sb, [j=J], sum(d[i,j]*y[i,j] for i in I) - y0[j] <= u*x[j])
        @constraint(sb, [i=I], sum(y[i,j] for j in J) == h_ξ[i,k])
    end

    SIPLIB.write_smps(model, filename)

    # Write .dro file

    # compute Wasserstein distance
    dist = zeros(nS,nS+nK)
    for s = S, ss = S
        dist[s,ss] = sqrt(norm(h[:,s] - h[:,ss])^2)
    end
    for s = S, k = K
        dist[s,nS+k] = sqrt(norm(h[:,s] - h_ξ[:,k])^2)
    end

    SIPLIB.write_wasserstein_dro(nS, nK, p, dist, 1.0, filename)

    return
end