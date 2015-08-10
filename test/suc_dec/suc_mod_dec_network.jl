# Julia script for unit commitment problem
# The original problem was written in AMPL scripts,
# which were provided by Changhyeok Lee
# Kibaek Kim - 2015 ANL MCS

# Based on suc_mod.jl: replaced load shift factor constraints by
# power flow and phase angle constraints with flow and angle.
# CT

# NOTE:
#   - Construct model in a local scope and feed it to DSP solver.
#   - This allows garbage collection to free memory.

include("clustering.jl");

# ---------------------------------
# Parameter setting at global scope
# ---------------------------------

if isdefined(:nScenarios) == false
	nScenarios  = 3;
end

# IEEE 118
nBuses      = 118;
nBranches   = 186;
nGenerators = 54;
nWinds      = 3;
# nPeriods    = 24;
# nPeriods    = 12;
nPeriods    = 6;
# nPeriods    = 1;
nSegments   = 4;
nCutEdges   = 0;

let

data_dir = "IEEE118";

# Relax integrality of variable Use
relax_integrality = false;
# relax_integrality = true;

# Apply graph partitioning to decompose problem
decomposition = true;
# decomposition = false;

# ---------------
# Read data files
# ---------------

tmpdat1 = readdlm(joinpath(data_dir, "generator.dat"), '\t');
tmpdat2 = readdlm(joinpath(data_dir, "generator_cost_function.dat"), '\t');
# tmpdat3 = readdlm(joinpath(data_dir, "load_profile_$(nPeriods).dat"), '\t');
tmpdat3 = readdlm(joinpath(data_dir, "load_profile_6_network.dat"), '\t');
tmpdat4 = readdlm(joinpath(data_dir, "load_distribution.dat"), '\t');
if isfile(joinpath(data_dir, "wind_profile.dat")) && isfile(joinpath(data_dir, "wind_distribution.dat"))
	tmpdat5 = readdlm(joinpath(data_dir, "wind_profile.dat"), '\t');
	tmpdat6 = readdlm(joinpath(data_dir, "wind_distribution.dat"), '\t');
end
tmpdat8 = readdlm(joinpath(data_dir, "branch.dat"), '\t');

# -----------------
# Parameter setting
# -----------------

BUSES      = 1:nBuses;
BRANCHES   = 1:nBranches;
GENERATORS = 1:nGenerators;
WINDS      = 1:nWinds;
PERIODS    = 1:nPeriods;
SEGMENTS   = 1:nSegments;
CUTEDGES   = 1:nCutEdges;

prob = ones(nScenarios) / nScenarios; # probabilities

#tmpdat1[:,1];                  # Index
gen_bus_id     = tmpdat1[:,2];  # Bus ID, where generators are attached
min_gen        = tmpdat1[:,3];  # Minimum power generation (MW)
max_gen        = tmpdat1[:,4];  # Maximum power generation (MW)
gen_0          = tmpdat1[:,5];  # Power generation in hour 0 (MW)
use_history    = tmpdat1[:,6];  # number of hours the generator has been on or off before hour 1
uptime         = tmpdat1[:,7];  # Minimum uptime (hours)
downtime       = tmpdat1[:,8];  # Minimum downtime (hours)
ramp_rate      = tmpdat1[:,9];  # Ramp rate per hour (MW/hour)
quick_start    = tmpdat1[:,12]; # Quick start capability?
fixed_cost_gen = tmpdat1[:,13]; # Fixed cost of running the generator for an hour ($)
cost_start     = tmpdat1[:,14]; # Cost of starting a generator ($)

# define quick starters
FASTGENS = sort(unique(GENERATORS .* int(quick_start)));
SLOWGENS = sort(unique(GENERATORS .* int(1 - quick_start)));
if FASTGENS[1] == 0
        FASTGENS = FASTGENS[2:length(FASTGENS)];
end
if SLOWGENS[1] == 0
        SLOWGENS = SLOWGENS[2:length(SLOWGENS)];
end

cost_gen     = tmpdat2[:,1:4]; # Marginal cost of production ($/MWh)
max_gen_sgmt = tmpdat2[:,5:8]; # Length of each power output segment (MW)

spin_resv_rate = 0.04; # Spinning reserve percentage
spin_notice    = 10;   # Spinning notice window (min)

total_demand = tmpdat3[:,1]; # Power load for each hour (MW)
demand_dist = zeros(nBuses); # Demand distribution
for i in 1:size(tmpdat4,1)
	demand_dist[tmpdat4[i,2]] = tmpdat4[i,3];
end

if nWinds > 0
	total_wind  = tmpdat5[:,1]; # Total wind power generation for each hour (MW)
	wind_bus_id = tmpdat6[:,2]; # Wind power generator ID
	wind_dist   = tmpdat6[:,3]; # Wind power distribution
else
	# Wind power is not used
	total_wind  = [];
	wind_bus_id = [];
	wind_dist   = [];
end

flowin_bus_id     = tmpdat8[:,2];            # flow-in bus ID (from)
flowout_bus_id    = tmpdat8[:,3];            # flow-out bus ID (to)
resistance        = tmpdat8[:,5];            # resistance
reactance         = tmpdat8[:,6];            # reactance
flow_max          = tmpdat8[:,8];            # Transmission line capacity

flowin_bus_id  = convert(Array{Int,1}, flowin_bus_id);
flowout_bus_id = convert(Array{Int,1}, flowout_bus_id);

# Adjust downtime/uptime according to number of periods
for i in 1:length(uptime)
	uptime[i] = ceil(uptime[i] * (nPeriods / 24.0))
end
for i in 1:length(downtime)
	downtime[i] = ceil(downtime[i] * (nPeriods / 24.0))
end

println("Done reading model file");

# -------------------------------
# Partition graph
# -------------------------------

if decomposition
	weights_type = :topology;
	# weights_type = :admittance;
	precision = 1e6; # METIS only accepts integer weights; if weight is a floating point, we multiply it by precision and round it

	weights = Dict{(Int,Int), Int}();
	if weights_type == :topology
		for (i,j) in zip(flowin_bus_id, flowout_bus_id)
			weights[(i,j)] = 1;
		end
	elseif weights_type == :admittance
		for (i,j,resist,react) in zip(flowin_bus_id, flowout_bus_id, resistance, reactance)
			float_weight = 1 / sqrt(resist^2 + react^2);
			weights[(i,j)] = round(Int, float_weight * precision);
		end
	end

	println("Partitioning graph...");
	@time graph, edgelist, partition = build_and_cluster_graph(nBuses, flowin_bus_id, flowout_bus_id, weights, nclusters);
	println("Done partitioning graph");

	# nclusters = maximum(partition);
	cut_edge_indices = collect_cut_edges(nBuses, edgelist, partition);

	# Create dummy nodes for decomposition

	nCutEdges = length(cut_edge_indices);
	println("$nCutEdges cut edges");

	# Dummy nodes correspond to cut edges (there may be more than one dummy node per node)
	# get_dummy_node_pair returns the pair of dummy node indices corresponding to a cut edge index
	originalNBuses = nBuses;
	get_dummy_nodes = function(i) # i is an index of cut_edge_indices
		((originalNBuses + 1) + 2 * (i - 1),
		 (originalNBuses + 1) + 2 * (i - 1) + 1)
	end

	# Update properties of dummy nodes
	nBuses += 2 * nCutEdges;
	for edge_index in cut_edge_indices
		edge = edgelist[edge_index];
		# This suffices to consider demand as demand is random on periods, not individual buses
		# push!(demand_dist, demand_dist[edge[1]]);
		# push!(demand_dist, demand_dist[edge[2]]);
		push!(demand_dist, 0);
		push!(demand_dist, 0);

		# Update partition to classify dummy nodes (useful for decomposition)
		push!(partition, partition[edge[1]]);
		push!(partition, partition[edge[2]]);
	end

	# We add two dummy edges per cut edge: we replace the original one and add a new one
	originalNBranches = nBranches;
	get_dummy_edges = function(i) # i is an index of cut_edge_indices
		(cut_edge_indices[i], originalNBranches + i)
	end

	# Update properties of dummy edges
	nBranches += nCutEdges;
	append!(flowin_bus_id, zeros(Int, nCutEdges));
	append!(flowout_bus_id, zeros(Int, nCutEdges));
	append!(resistance, zeros(Int, nCutEdges));
	append!(reactance, zeros(Int, nCutEdges));
	append!(flow_max, zeros(Int, nCutEdges));
	for (i, edge_index) in enumerate(cut_edge_indices)
		edge = edgelist[edge_index];
		dummy_node_pair = get_dummy_nodes(i);
		dummy_edge_pair = get_dummy_edges(i);

		flowin_bus_id[dummy_edge_pair[1]] = edge[1];
		flowout_bus_id[dummy_edge_pair[1]] = dummy_node_pair[1];
		flowin_bus_id[dummy_edge_pair[2]] = dummy_node_pair[2];
		flowout_bus_id[dummy_edge_pair[2]] = edge[2];

		resistance[dummy_edge_pair[2]] = resistance[dummy_edge_pair[1]];
		reactance[dummy_edge_pair[2]] = reactance[dummy_edge_pair[1]];
		flow_max[dummy_edge_pair[2]] = flow_max[dummy_edge_pair[1]];
	end


	ORIGINALBRANCHES = 1:originalNBranches;
	ORIGINALBUSES    = 1:originalNBuses;
	BUSES            = 1:nBuses;
	BRANCHES         = 1:nBranches;
	CUTEDGES         = 1:nCutEdges;

	# Plot modified graph
	# graph, edgelist = graph_from_data(nBuses, flowin_bus_id, flowout_bus_id);
	# plot_graph(nBuses, edgelist);
end

# -------------------------------
# Initialize auxiliary parameters
# -------------------------------

println("Constructing model...");

use_0         = zeros(nGenerators);      # Unit commitment in hour 0
downtime_init = zeros(nGenerators);      # Initial minimum downtime
uptime_init   = zeros(nGenerators);      # Initial minimum uptime
demand        = zeros(nBuses, nPeriods); # Power load
wind          = zeros(nWinds, nPeriods); # Wind power generation

for i in GENERATORS
	if use_history[i] > 0
		use_0[i] = 1;
	else
		use_0[i] = 0;
	end
	downtime_init[i] = max(0, downtime[i] + use_history[i]) * (1 - use_0[i]);
	uptime_init[i]   = max(0, uptime[i] - use_history[i]) * use_0[i];
end

for t in PERIODS
	for n in BUSES
		demand[n,t] = demand_dist[n] / sum(demand_dist) * total_demand[t];
	end
	for n in WINDS
		wind[n,t] = wind_dist[n] / sum(wind_dist) * total_wind[t];
	end
end

# Susceptance of the lines
susceptance = -reactance ./ (resistance.^2 + reactance.^2);

# -------------------
# Scenario generation
# -------------------

srand(1);
load_uncertainty = 0.03; # Load uncertainty
wind_uncertainty = 0.55;  # Wind power generation uncertainty

# Load scenarios
total_demand_scen = zeros(nPeriods, nScenarios);
demand_scen = zeros(nBuses, nPeriods, nScenarios);
total_demand_scen[:,1] = total_demand;
#demand_scen[:,:,1]     = demand;
for s in 1:nScenarios
	total_demand_scen[:,s] = total_demand * (1 - load_uncertainty) + rand(nPeriods) .* total_demand * load_uncertainty * 2;
	demand_scen[:,:,s] = demand_dist ./ sum(demand_dist) * total_demand_scen[:,s]';
end

# Wind power scenarios
wind_scen = zeros(nWinds, nPeriods, nScenarios);
#wind_scen[:,:,1] = wind;
for s in 1:nScenarios
	wind_scen[:,:,s] = wind * (1 - wind_uncertainty) + rand(nWinds, nPeriods) .* wind * wind_uncertainty * 2;
end
total_wind_scen = reshape(sum(wind_scen,1), nPeriods, nScenarios);



# -----------------
# Release file data
# -----------------

tmpdat1 = 0;
tmpdat2 = 0;
tmpdat3 = 0;
tmpdat4 = 0;
tmpdat5 = 0;
tmpdat6 = 0;
tmpdat8 = 0;

# ----------------
# StochJuMP object
# ----------------
m = StochasticModel(nScenarios);

# ------------------------------------------------------
# The following parameters need to be defined with data.
# ------------------------------------------------------
# BUSES
# BRANCHES
# GENERATORS
# WINDS
# PERIODS
# SEGMENTS
# gen_bus_id       [GENERATORS]
# cost_start       [GENERATORS]
# fixed_cost_gen   [GENERATORS]
# cost_gen         [GENERATORS, SEGMENTS]
# use_history      [GENERATORS]
# downtime         [GENERATORS]
# uptime           [GENERATORS]
# use_0            [GENERATORS]
# downtime_init    [GENERATORS]
# uptime_init      [GENERATORS]
# min_gen          [GENERATORS]
# max_gen          [GENERATORS]
# max_gen_sgmt     [GENERATORS, SEGMENTS]
# ramp_rate        [GENERATORS]
# gen_0            [GENERATORS]
# spin_resv_rate
# spin_notice
# total_demand     [PERIODS]
# demand           [BUSES,PERIODS]
# wind_bus_id      [WINDS]
# total_wind       [PERIODS]
# wind             [WINDS,PERIODS]
# resistance       [BRANCHES]
# reactance        [BRANCHES]
# flow_max         [BRANCHES]

# ---------------------
# First-stage Variables
# ---------------------
if relax_integrality
	@defVar(m, 0 <= Use[i=SLOWGENS, t=PERIODS] <= 1)  # Generator on/off indicator
else
	@defVar(m, Use[i=SLOWGENS, t=PERIODS], Bin)       # Generator on/off indicator
end
@defVar(m, 0 <= Up[i=SLOWGENS, t=PERIODS] <= 1)   # Start up indicator
@defVar(m, 0 <= Down[i=SLOWGENS, t=PERIODS] <= 1) # Shut down indicator

# ------------------------------
# First-stage Objective function
# ------------------------------
@setObjective(m, Min,
	sum{cost_start[i] * Up[i,t], i=SLOWGENS, t=PERIODS}
	+ sum{fixed_cost_gen[i] * Use[i,t], i=SLOWGENS, t=PERIODS})

# -----------------------
# First-stage Constraints
# -----------------------

# Linking Use / Up / Down variables
@addConstraint(m, LINKING_SHUT_DOWN0[i=SLOWGENS],
	Down[i,1] <= use_0[i])
@addConstraint(m, LINKING_SHUT_DOWN[i=SLOWGENS, t=2:nPeriods],
	Use[i,t-1] >= Down[i,t])
@addConstraint(m, LINKING_START_UP0[i=SLOWGENS],
	Up[i,1] <= 1 - use_0[i])
@addConstraint(m, LINKING_START_UP[i=SLOWGENS, t=2:nPeriods],
	1 - Use[i,t-1] >= Up[i,t])
@addConstraint(m, LINKING_BOTH0[i=SLOWGENS],
	Use[i,1] - use_0[i] == Up[i,1] - Down[i,1])
@addConstraint(m, LINKING_BOTH[i=SLOWGENS, t=2:nPeriods],
	Use[i,t] - Use[i,t-1] == Up[i,t] - Down[i,t])

# Min down time
@addConstraint(m, MIN_DOWN_INIT[i=SLOWGENS, t=1:min(downtime_init[i],nPeriods)],
	Use[i,t] == 0)
@addConstraint(m, MIN_DOWN_S1[i=SLOWGENS, t=PERIODS, s=max(1,t-downtime[i]+1):t],
	1 - Use[i,t] >= Down[i,s])
@addConstraint(m, MIN_DOWN_S2[i=SLOWGENS, t=PERIODS],
	1 - Use[i,t] >= sum{Down[i,s], s=max(1,t-downtime[i]+1):t})

# Min up time
@addConstraint(m, MIN_UP_INIT[i=SLOWGENS, t=1:min(uptime_init[i],nPeriods)],
	Use[i,t] == 1)
@addConstraint(m, MIN_UP_S1[i=SLOWGENS, t=PERIODS, s=max(1,t-uptime[i]+1):t],
	Use[i,t] >= Up[i,s])
@addConstraint(m, MIN_UP_S2[i=SLOWGENS, t=PERIODS],
	Use[i,t] >= sum{Up[i,s], s=max(1,t-uptime[i]+1):t})

# -----------------
# For each scenario
# -----------------
@second_stage m s begin

	# ----------------
	# Stochastic block
	# ----------------
	sb = StochasticBlock(m, prob[s])

	# ----------------------
	# Second-stage Variables
	# ----------------------

	if relax_integrality
	    @defVar(sb, 0 <= UseF[i=FASTGENS, t=PERIODS] <= 1)        # Generator on/off indicator
	else
	    @defVar(sb, UseF[i=FASTGENS, t=PERIODS], Bin)        # Generator on/off indicator
	end
    @defVar(sb, 0 <= UpF[i=FASTGENS, t=PERIODS] <= 1)    # Start up indicator
    @defVar(sb, 0 <= DownF[i=FASTGENS, t=PERIODS] <= 1)  # Shut down indicator

	@defVar(sb, 0 <= Gen[i=GENERATORS, t=PERIODS] <= max_gen[i])       # Power generation
	@defVar(sb, 0 <= Gen_Sgmt[i=GENERATORS, k=SEGMENTS, t=PERIODS] <= max_gen_sgmt[i,k])
	@defVar(sb, 0 <= Spin_Resv[i=GENERATORS, t=PERIODS] <= spin_notice / 60. * ramp_rate[i]) # Spinning reserve
	@defVar(sb, -flow_max[l] <= Flow[l=BRANCHES, t=PERIODS] <= flow_max[l]);              # power flow (Transmission line constraints)
	@defVar(sb, -360.0 <= Phase[n=BUSES, t=PERIODS] <= 360.0);                            # Phase angle

	# -------------------------------
	# Second-stage Objective function
	# -------------------------------

	@setObjective(sb, Min,
        sum{cost_start[i] * UpF[i,t], i=FASTGENS, t=PERIODS}
        + sum{fixed_cost_gen[i] * UseF[i,t], i=FASTGENS, t=PERIODS}
		+ sum{cost_gen[i,k] * Gen_Sgmt[i,k,t], i=GENERATORS, k=SEGMENTS, t=PERIODS})

	# ------------------------
	# Second-stage Constraints
	# ------------------------

    # Linking Use / Up / Down variables
    @addConstraint(sb, FAST_LINKING_SHUT_DOWN0[i=FASTGENS],
            DownF[i,1] <= use_0[i])
    @addConstraint(sb, FAST_LINKING_SHUT_DOWN[i=FASTGENS, t=2:nPeriods],
            UseF[i,t-1] >= DownF[i,t])
    @addConstraint(sb, FAST_LINKING_START_UP0[i=FASTGENS],
            UpF[i,1] <= 1 - use_0[i])
    @addConstraint(sb, FAST_LINKING_START_UP[i=FASTGENS, t=2:nPeriods],
            1 - UseF[i,t-1] >= UpF[i,t])
    @addConstraint(sb, FAST_LINKING_BOTH0[i=FASTGENS],
            UseF[i,1] - use_0[i] == UpF[i,1] - DownF[i,1])
    @addConstraint(sb, FAST_LINKING_BOTH[i=FASTGENS, t=2:nPeriods],
            UseF[i,t] - UseF[i,t-1] == UpF[i,t] - DownF[i,t])

    # Min down time
    @addConstraint(sb, FAST_MIN_DOWN_INIT[i=FASTGENS, t=1:min(downtime_init[i],nPeriods)],
            UseF[i,t] == 0)
    @addConstraint(sb, FAST_MIN_DOWN_S1[i=FASTGENS, t=PERIODS, tt=max(1,t-downtime[i]+1):t],
            1 - UseF[i,t] >= DownF[i,tt])
    @addConstraint(sb, FAST_MIN_DOWN_S2[i=FASTGENS, t=PERIODS],
            1 - UseF[i,t] >= sum{DownF[i,tt], tt=max(1,t-downtime[i]+1):t})

    # Min up time
    @addConstraint(sb, FAST_MIN_UP_INIT[i=FASTGENS, t=1:min(uptime_init[i],nPeriods)],
            UseF[i,t] == 1)
    @addConstraint(sb, FAST_MIN_UP_S1[i=FASTGENS, t=PERIODS, tt=max(1,t-uptime[i]+1):t],
            UseF[i,t] >= UpF[i,tt])
    @addConstraint(sb, FAST_MIN_UP_S2[i=FASTGENS, t=PERIODS],
            UseF[i,t] >= sum{UpF[i,tt], tt=max(1,t-uptime[i]+1):t})

	# Ramping rate in normal operating status
	@addConstraint(sb, RAMP_DOWN0[i=GENERATORS],
		gen_0[i] - Gen[i,1] <= ramp_rate[i])
	@addConstraint(sb, RAMP_DOWN[i=GENERATORS, t=2:nPeriods],
		Gen[i,t-1] - Gen[i,t] <= ramp_rate[i])
	@addConstraint(sb, RAMP_UP0[i=GENERATORS],
		Gen[i,1] - gen_0[i] + Spin_Resv[i,1] <= ramp_rate[i])
	@addConstraint(sb, RAMP_UP[i=GENERATORS, t=2:nPeriods],
		Gen[i,t] - Gen[i,t-1] + Spin_Resv[i,t] <= ramp_rate[i])

	# Spinning reserve requirement for system
	# **** TODO **** This is a global constraint that we are currently ignoring. Might be ok to split it up?
	# @addConstraint(sb, SPIN_RESV_REQ[t=PERIODS],
	# 	sum{Spin_Resv[i,t], i=GENERATORS}
	# 	>= spin_resv_rate * (total_demand_scen[t,s] - total_wind_scen[t,s]))

	# Spinning reserve capacity for individual unit
    @addConstraint(sb, SPIN_RESV_MAX_SLOW[i=SLOWGENS, t=PERIODS],
            Spin_Resv[i,t] <= spin_notice / 60. * ramp_rate[i] * Use[i,t])
    @addConstraint(sb, SPIN_RESV_MAX_FAST[i=FASTGENS, t=PERIODS],
            Spin_Resv[i,t] <= spin_notice / 60. * ramp_rate[i] * UseF[i,t])

	# Power output capacity constraints
    @addConstraint(sb, POWER_OUTPUT_SLOW[i=SLOWGENS, t=PERIODS],
            Gen[i,t] == min_gen[i] * Use[i,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
    @addConstraint(sb, POWER_SEGMENT_SLOW[i=SLOWGENS, k=SEGMENTS, t=PERIODS],
            Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * Use[i,t])
    @addConstraint(sb, POWER_MAX_SLOW[i=SLOWGENS, t=PERIODS],
            Gen[i,t] + Spin_Resv[i,t] <= max_gen[i] * Use[i,t])
    @addConstraint(sb, POWER_OUTPUT_FAST[i=FASTGENS, t=PERIODS],
            Gen[i,t] == min_gen[i] * UseF[i,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
    @addConstraint(sb, POWER_SEGMENT_FAST[i=FASTGENS, k=SEGMENTS, t=PERIODS],
            Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * UseF[i,t])
    @addConstraint(sb, POWER_MAX_FAST[i=FASTGENS, t=PERIODS],
            Gen[i,t] + Spin_Resv[i,t] <= max_gen[i] * UseF[i,t])

	if decomposition
		# Power balance constraints for system
		@addConstraint(sb, POWER_BALANCE[n=ORIGINALBUSES, t=PERIODS],
			sum{Flow[l,t], l=BRANCHES; flowout_bus_id[l] == n}   # incoming arcs (destination is n)
			- sum{Flow[l,t], l=BRANCHES; flowin_bus_id[l] == n}  # outgoing arcs (source is n)
			+ sum{Gen[i,t], i=GENERATORS; gen_bus_id[i] == n}
			== demand_scen[n,t,s] - sum{wind_scen[i,t,s], i=WINDS; wind_bus_id[i] == n});
	else
		# Power balance constraints for system
		@addConstraint(sb, POWER_BALANCE[n=BUSES, t=PERIODS],
			sum{Flow[l,t], l=BRANCHES; flowout_bus_id[l] == n}   # incoming arcs (destination is n)
			- sum{Flow[l,t], l=BRANCHES; flowin_bus_id[l] == n}  # outgoing arcs (source is n)
			+ sum{Gen[i,t], i=GENERATORS; gen_bus_id[i] == n}
			== demand_scen[n,t,s] - sum{wind_scen[i,t,s], i=WINDS; wind_bus_id[i] == n});
	end

	# Power flow equation
	@addConstraint(sb, POWER_FLOW[l=BRANCHES, t=PERIODS],
		Flow[l,t] == susceptance[l] * (Phase[flowout_bus_id[l],t] - Phase[flowin_bus_id[l],t]));

	# Angle-based decomposition
	# if decomposition
	# 	# Connects different clusters through phase angle
	# 	@addConstraint(sb, INTERCLASS_PHASE_1[l=CUTEDGES, t=PERIODS],
	# 		Phase[get_dummy_nodes(l)[1],t] == Phase[edgelist[cut_edge_indices[l]][2],t]);
	# 	@addConstraint(sb, INTERCLASS_PHASE_2[l=CUTEDGES, t=PERIODS],
	# 		Phase[get_dummy_nodes(l)[2],t] == Phase[edgelist[cut_edge_indices[l]][1],t]);
	# end

	# Angle-based decomposition
	if decomposition
		for l in CUTEDGES
			for t in PERIODS
				DSPsolver.addCouplingConstraint(sb, @JuMP.LinearConstraint(
					Phase[get_dummy_nodes(l)[1],t] == Phase[edgelist[cut_edge_indices[l]][2],t]));
				DSPsolver.addCouplingConstraint(sb, @JuMP.LinearConstraint(
					Phase[get_dummy_nodes(l)[2],t] == Phase[edgelist[cut_edge_indices[l]][1],t]));
			end
		end
	end

end

if decomposition
	for i in SLOWGENS, t in PERIODS
		DSPsolver.setVarSubproblem(m, Use[i,t], partition[gen_bus_id[i]]);
		DSPsolver.setVarSubproblem(m, Up[i,t], partition[gen_bus_id[i]]);
		DSPsolver.setVarSubproblem(m, Down[i,t], partition[gen_bus_id[i]]);
	end
	for sb in getchildren(m)
		UseF = JuMP.getVar(sb, :UseF);
		UpF = JuMP.getVar(sb, :UpF);
		DownF = JuMP.getVar(sb, :DownF);
		Gen = JuMP.getVar(sb, :Gen);
		Gen_Sgmt = JuMP.getVar(sb, :Gen_Sgmt);
		Spin_Resv = JuMP.getVar(sb, :Spin_Resv);
		Flow = JuMP.getVar(sb, :Flow);
		Phase = JuMP.getVar(sb, :Phase);

		for i in FASTGENS, t in PERIODS
			DSPsolver.setVarSubproblem(sb, UseF[i,t], partition[gen_bus_id[i]]);
			DSPsolver.setVarSubproblem(sb, UpF[i,t], partition[gen_bus_id[i]]);
			DSPsolver.setVarSubproblem(sb, DownF[i,t], partition[gen_bus_id[i]]);
		end
		for i in GENERATORS, t in PERIODS
			DSPsolver.setVarSubproblem(sb, Gen[i,t], partition[gen_bus_id[i]]);
			DSPsolver.setVarSubproblem(sb, Spin_Resv[i,t], partition[gen_bus_id[i]]);
		end
		for i in GENERATORS, k in SEGMENTS, t in PERIODS
			DSPsolver.setVarSubproblem(sb, Gen_Sgmt[i,k,t], partition[gen_bus_id[i]]);
		end
		for l in BRANCHES, t in PERIODS
			if partition[flowout_bus_id[l]] != partition[flowin_bus_id[l]]
				println("Error: Network not decomposed");
				exit();
			end
			DSPsolver.setVarSubproblem(sb, Flow[l,t], partition[flowout_bus_id[l]]);
		end
		for n in BUSES, t in PERIODS
			DSPsolver.setVarSubproblem(sb, Phase[n,t], partition[n]);
		end
	end
end

# Load data to DSP
DSPsolver.loadProblem(m);

println("Done constructing model");

end # End of let
