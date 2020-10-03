using SIPLIB
using StructJuMP, JuMP
using LinearAlgebra
using Random

DRDCAP(nR::Int, nN::Int, nT::Int, nS::Int, nK::Int, ϵ::Float64, filename::String) = DRDCAP(nR, nN, nT, nS, nK, 0, ϵ, filename)

function DRDCAP(nR::Int, nN::Int, nT::Int, nS::Int, nK::Int, nE::Int, ϵ::Float64, filename::String)

    # set random seed (default=1)
    Random.seed!(1)

    # generate & store instance data
    ## sets
    R = 1:nR
    N = 1:nN
    T = 1:nT
    S = 1:nS # number of references
    K = 1:nK # number of discretization points

    ## parameters
    a = rand(nR, nT) * 5 .+ 5
    b = rand(nR, nT) * 40 .+ 10
    c = rand(nR, nN, nT, nS) * 2 .+ 5  # (5,7]
    c0 = rand(nN, nT, nS) * 200 .+ 500 # (500,700]
    d = rand(nN, nT, nS) * 0.5 .+ 1  # (1.0,1.5]
    p = ones(nS) / nS
    Pr = 1.0 / (nS+nK)

    # generate many scenarios ahead
    # This allows us to consistently increase the number of scenarios for experiment.
    max_nK = Int(1e+6)
    c_ξ = rand(nR, nN, nT, max_nK) * 5 .+ 5  # (5,10]
    c0_ξ = rand(nN, nT, max_nK) * 500 .+ 500 # (500,1000]
    d_ξ = rand(nN, nT, max_nK) .+ 0.5        # (0.5,1.5]

    # construct JuMP.Model
    model = StructuredModel(num_scenarios = nS+nK)

    ## 1st stage
    if nE > 0
        @variable(model, x2[i=R,t=T,e=1:nE], Bin)
        @expression(model, x[i=R,t=T], sum(x2[i,t,e]/2^(e-1) for e=1:nE))
    else
        @variable(model, x[i=R,t=T] >= 0)
    end
    @variable(model, u[i=R,t=T], Bin)
    @objective(model, Min, sum(a[i,t]*x[i,t] + b[i,t]*u[i,t] for i in R for t in T))
    @constraint(model, [i=R,t=T], x[i,t] - u[i,t] <= 0)

    ## 2nd stage
    for ss = S
        sb = StructuredModel(parent=model, id = ss, prob = Pr)
        @variable(sb, y[i=0:nR, j=N, t=T], Bin)
        @objective(sb, Min, 
              sum(c0[j,t,ss] * y[0,j,t] for j=N, t=T) 
            + sum(c[i,j,t,ss] * y[i,j,t] for i=R, j=N, t=T))
        @constraint(sb, [i=R,t=T], -sum(x[i,tau] for tau in 1:t) + sum(d[j,t,ss]*y[i,j,t] for j in N) <= 0)
        @constraint(sb, [j=N,t=T], sum(y[i,j,t] for i in 0:nR) == 1)
    end
    for k = K
        sb = StructuredModel(parent=model, id = nS+k, prob = Pr)
        @variable(sb, y[i=0:nR, j=N, t=T], Bin)
        @objective(sb, Min, 
              sum(c0_ξ[j,t,k] * y[0,j,t] for j=N, t=T) 
            + sum(c_ξ[i,j,t,k] * y[i,j,t] for i=R, j=N, t=T))
        @constraint(sb, [i=R,t=T], -sum(x[i,tau] for tau in 1:t) + sum(d_ξ[j,t,k]*y[i,j,t] for j in N) <= 0)
        @constraint(sb, [j=N,t=T], sum(y[i,j,t] for i in 0:nR) == 1)
    end

    SIPLIB.write_smps(model, filename)

    # Write .dro file
    dist = zeros(nS,nS+nK)
    for s = S, ss = S
        dist[s,ss] = sqrt(norm(c[:,:,:,s] - c[:,:,:,ss])^2 + norm(c0[:,:,s] - c0[:,:,ss])^2 + norm(d[:,:,s] - d[:,:,ss])^2)
    end
    for s = S, k = K
        dist[s,nS+k] = sqrt(norm(c[:,:,:,s] - c_ξ[:,:,:,k])^2 + norm(c0[:,:,s] - c0_ξ[:,:,k])^2 + norm(d[:,:,s] - d_ξ[:,:,k])^2)
    end

    SIPLIB.write_wasserstein_dro(nS, nK, p, dist, ϵ, filename)

	"""
	Deterministic equivalent model
	"""
#=
    # construct JuMP.Model
    detm = Model()

    ## 1st stage
    @variable(detm, x_d[i=R,t=T] >= 0)
    @variable(detm, u_d[i=R,t=T], Bin)
    @variable(detm, α_d >= 0)
    @variable(detm, β_d[s=S])
    @variable(detm, y_d[i=0:nR, j=N, t=T, k=1:(nS+nK)], Bin)
    @objective(detm, Min, 
        sum(a[i,t]*x_d[i,t] + b[i,t]*u_d[i,t] for i in R for t in T)
        + ϵ * α_d 
        + sum(p[s] * β_d[s] for s=S))
    @constraint(detm, [i=R,t=T], x_d[i,t] - u_d[i,t] <= 0)

    ## 2nd stage
    for ss = S
        @constraint(detm, [s=S], 
            dist[s,ss] * α_d + β_d[s] - sum(c0[j,t,ss] * y_d[0,j,t,ss] for j=N, t=T) - sum(c[i,j,t,ss] * y_d[i,j,t,ss] for i=R, j=N, t=T) >= 0)
        @constraint(detm, [i=R,t=T], -sum(x_d[i,tau] for tau in 1:t) + sum(d[j,t,ss]*y_d[i,j,t,ss] for j in N) <= 0)
        @constraint(detm, [j=N,t=T], sum(y_d[i,j,t,ss] for i in 0:nR) == 1)
    end
    for k = K
        @constraint(detm, [s=S], 
            dist[s,nS+k] * α_d + β_d[s] - sum(c0_ξ[j,t,k] * y_d[0,j,t,nS+k] for j=N, t=T) - sum(c_ξ[i,j,t,k] * y_d[i,j,t,nS+k] for i=R, j=N, t=T) >= 0)
        @constraint(detm, [i=R,t=T], -sum(x_d[i,tau] for tau in 1:t) + sum(d_ξ[j,t,k]*y_d[i,j,t,nS+k] for j in N) <= 0)
        @constraint(detm, [j=N,t=T], sum(y_d[i,j,t,nS+k] for i in 0:nR) == 1)
    end

	mof_model = MathOptFormat.MPS.Model()
	MOI.copy_to(mof_model, backend(detm))
	MOI.write_to_file(mof_model, "$filename.mps")
=#

    return
end

function DRDCAP_as_SMIP(nR::Int, nN::Int, nT::Int, nS::Int, nK::Int, ϵ::Float64, filename::String)

    # set random seed (default=1)
    Random.seed!(1)

    # generate & store instance data
    ## sets
    R = 1:nR
    N = 1:nN
    T = 1:nT
    S = 1:nS
    K = 1:nK

    ## parameters
    a = rand(nR, nT) * 5 .+ 5
    b = rand(nR, nT) * 40 .+ 10
    c = rand(nR, nN, nT, nS) * 2 .+ 5  # (5,7]
    c0 = rand(nN, nT, nS) * 200 .+ 500 # (500,700]
    d = rand(nN, nT, nS) * 0.5 .+ 1  # (1.0,1.5]
    p = ones(nS) / nS
    Pr = 1.0 / (nS+nK)

    # generate many scenarios ahead
    # This allows us to consistently increase the number of scenarios for experiment.
    max_nK = Int(1e+6)
    c_ξ = rand(nR, nN, nT, max_nK) * 5 .+ 5  # (5,10]
    c0_ξ = rand(nN, nT, max_nK) * 500 .+ 500 # (500,1000]
    d_ξ = rand(nN, nT, max_nK) .+ 0.5        # (0.5,1.5]

    dist = zeros(nS,nS+nK)
    for s = S, ss = S
        dist[s,ss] = sqrt(norm(c[:,:,:,s] - c[:,:,:,ss])^2 + norm(c0[:,:,s] - c0[:,:,ss])^2 + norm(d[:,:,s] - d[:,:,ss])^2)
    end
    for s = S, k = K
        dist[s,nS+k] = sqrt(norm(c[:,:,:,s] - c_ξ[:,:,:,k])^2 + norm(c0[:,:,s] - c0_ξ[:,:,k])^2 + norm(d[:,:,s] - d_ξ[:,:,k])^2)
    end

    # construct JuMP.Model
    model = StructuredModel(num_scenarios = nS+nK)

    ## 1st stage
    @variable(model, x[i=R,t=T] >= 0)
    @variable(model, u[i=R,t=T], Bin)
    @variable(model, 0 <= α <= 1e+6)
    @variable(model, -1e+6 <= β[s=S] <= 1e+6)
    @objective(model, Min, 
        sum(a[i,t]*x[i,t] + b[i,t]*u[i,t] for i in R for t in T)
        + ϵ * α 
        + sum(p[s] * β[s] for s=S))
    @constraint(model, [i=R,t=T], x[i,t] - u[i,t] <= 0)

    ## 2nd stage
    for ss = S
        sb = StructuredModel(parent=model, id = ss, prob = Pr)
        @variable(sb, y[i=0:nR, j=N, t=T], Bin)
        @constraint(sb, [s=S], 
            dist[s,ss] * α + β[s] - sum(c0[j,t,ss] * y[0,j,t] for j=N, t=T) - sum(c[i,j,t,ss] * y[i,j,t] for i=R, j=N, t=T) >= 0)
        @constraint(sb, [i=R,t=T], -sum(x[i,tau] for tau in 1:t) + sum(d[j,t,ss]*y[i,j,t] for j in N) <= 0)
        @constraint(sb, [j=N,t=T], sum(y[i,j,t] for i in 0:nR) == 1)
    end
    for k = K
        sb = StructuredModel(parent=model, id = nS+k, prob = Pr)
        @variable(sb, y[i=0:nR, j=N, t=T], Bin)
        @constraint(sb, [s=S], 
            dist[s,nS+k] * α + β[s] - sum(c0_ξ[j,t,k] * y[0,j,t] for j=N, t=T) - sum(c_ξ[i,j,t,k] * y[i,j,t] for i=R, j=N, t=T) >= 0)
        @constraint(sb, [i=R,t=T], -sum(x[i,tau] for tau in 1:t) + sum(d_ξ[j,t,k]*y[i,j,t] for j in N) <= 0)
        @constraint(sb, [j=N,t=T], sum(y[i,j,t] for i in 0:nR) == 1)
    end

    SIPLIB.write_smps(model, filename)

    return
end

#=
if length(ARGS) != 5
	error("invalid number of arguments")
end

m = parse(Int, ARGS[1])
n = parse(Int, ARGS[2])
T = parse(Int, ARGS[3])
s = parse(Int, ARGS[4])
k = parse(Int, ARGS[5])
DRDCAP(m, n, T, s, k, 0.0, "drdcap_$(m)$(n)$(T)_$(s)_$(k)")
=#

