# Julia script for unit commitment problem
# - saving as SMPS files
# Kibaek Kim - 2016 ANL-MCS

using JuMP

nScenarios = parse(Int,ARGS[1]);

# ---------------
# Read data files
# ---------------

tmpdat1 = readdlm("IEEE118/generator.dat", '\t');
tmpdat2 = readdlm("IEEE118/generator_cost_function.dat", '\t');
tmpdat3 = readdlm("IEEE118/load_profile.dat", '\t');
tmpdat4 = readdlm("IEEE118/load_distribution.dat", '\t');
tmpdat6 = readdlm("IEEE118/wind_distribution.dat", '\t');
tmpdat7 = readdlm("IEEE118/shift_factor.dat", '\t');
tmpdat8 = readdlm("IEEE118/branch.dat", '\t');

# -----------------
# Parameter setting
# -----------------

if isdefined(:nScenarios) == false
        nScenarios  = 3;
end
nBuses      = 118;
nBranches   = 186;
nGenerators = 54;
nWinds      = 3;
nPeriods    = 24;
nSegments   = 4;

BUSES      = 1:nBuses;
BRANCHES   = 1:nBranches;
GENERATORS = 1:nGenerators;
WINDS      = 1:nWinds;
PERIODS    = 1:nPeriods;
SEGMENTS   = 1:nSegments;

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
FASTGENS = sort(unique(GENERATORS .* round(Int64,quick_start)));
SLOWGENS = sort(unique(GENERATORS .* round(Int64,1 - quick_start)));
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
        demand_dist[round(Int,tmpdat4[i,2])] = tmpdat4[i,3];
end

wind_bus_id = tmpdat6[:,2]; # Wind power generator ID

# Load shift factor (This is pre-computed using branch_from_bus and branch_to_bus.)
load_shift_factor = tmpdat7';
flow_max          = tmpdat8[:,8]; # Transmission line capacity

# -------------------------------
# Initialize auxiliary parameters
# -------------------------------

use_0         = zeros(nGenerators);      # Unit commitment in hour 0
downtime_init = zeros(nGenerators);      # Initial minimum downtime
uptime_init   = zeros(nGenerators);      # Initial minimum uptime
demand        = zeros(nBuses, nPeriods); # Power load

for i in GENERATORS
        if use_history[i] > 0
                use_0[i] = 1;
        else
                use_0[i] = 0;
        end
	downtime[i] = round(Int, downtime[i]);
	uptime[i]   = round(Int, uptime[i]);
        downtime_init[i] = round(Int, max(0, downtime[i] + use_history[i]) * (1 - use_0[i]));
        uptime_init[i]   = round(Int, max(0, uptime[i] - use_history[i]) * use_0[i]);
end

for t in PERIODS
        for n in BUSES
                demand[n,t] = demand_dist[n] / sum(demand_dist) * total_demand[t];
        end
end

# -------------------
# Scenario generation
# -------------------

# Wind power scenarios
wind_scen = zeros(nWinds, nPeriods, nScenarios);
for i in 1:nWinds
        winddat = readdlm(string("IEEE118/wp64-", i, ".txt"));
        wind_scen[i,:,:] = winddat[:,1:nScenarios];
end
total_wind_scen = reshape(sum(wind_scen,1), nPeriods, nScenarios);

# ----------------
# StochJuMP object
# ----------------
m = Model();

# ---------------------
# First-stage Variables
# ---------------------
@defVar(m, Use[i=SLOWGENS, t=PERIODS], Bin);         # Generator on/off indicator
@defVar(m, 0 <= Up[i=SLOWGENS, t=PERIODS] <= 1);     # Start up indicator
@defVar(m, 0 <= Down[i=SLOWGENS, t=PERIODS] <= 1);   # Shut down indicator

# ----------------------
# Second-stage Variables
# ----------------------

startCol4 = MathProgBase.numvar(m)+1;

@defVar(m, UseF[i=FASTGENS, t=PERIODS], Bin);        # Generator on/off indicator
@defVar(m, 0 <= UpF[i=FASTGENS, t=PERIODS] <= 1);    # Start up indicator
@defVar(m, 0 <= DownF[i=FASTGENS, t=PERIODS] <= 1);  # Shut down indicator
@defVar(m, 0 <= Gen[i=GENERATORS, t=PERIODS] <= max_gen[i]);                             # Power generation
@defVar(m, 0 <= Gen_Sgmt[i=GENERATORS, k=SEGMENTS, t=PERIODS] <= max_gen_sgmt[i,k]);     # generation cost sement 
@defVar(m, 0 <= Spin_Resv[i=GENERATORS, t=PERIODS] <= spin_notice / 60. * ramp_rate[i]); # Spinning reserve
@defVar(m, -flow_max[l] <= Flow[l=BRANCHES, t=PERIODS] <= flow_max[l]);                  # flow amount

# ------------------
# Objective function
# ------------------
@setObjective(m, Min,
        sum{cost_start[i] * Up[i,t], i=SLOWGENS, t=PERIODS}
        + sum{fixed_cost_gen[i] * Use[i,t], i=SLOWGENS, t=PERIODS}
        + sum{cost_start[i] * UpF[i,t], i=FASTGENS, t=PERIODS}
        + sum{fixed_cost_gen[i] * UseF[i,t], i=FASTGENS, t=PERIODS}
        + sum{cost_gen[i,k] * Gen_Sgmt[i,k,t], i=GENERATORS, k=SEGMENTS, t=PERIODS});

# -----------
# Constraints
# -----------

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


# ------------------------
# Second-stage constraints
# ------------------------

startRow4 = MathProgBase.numlinconstr(m)+1;

# Linking Use / Up / Down variables
@addConstraint(m, FAST_LINKING_SHUT_DOWN0[i=FASTGENS],
        DownF[i,1] <= use_0[i])
@addConstraint(m, FAST_LINKING_SHUT_DOWN[i=FASTGENS, t=2:nPeriods],
        UseF[i,t-1] >= DownF[i,t])
@addConstraint(m, FAST_LINKING_START_UP0[i=FASTGENS],
        UpF[i,1] <= 1 - use_0[i])
@addConstraint(m, FAST_LINKING_START_UP[i=FASTGENS, t=2:nPeriods],
        1 - UseF[i,t-1] >= UpF[i,t])
@addConstraint(m, FAST_LINKING_BOTH0[i=FASTGENS],
        UseF[i,1] - use_0[i] == UpF[i,1] - DownF[i,1])
@addConstraint(m, FAST_LINKING_BOTH[i=FASTGENS, t=2:nPeriods],
        UseF[i,t] - UseF[i,t-1] == UpF[i,t] - DownF[i,t])

# Min down time
@addConstraint(m, FAST_MIN_DOWN_INIT[i=FASTGENS, t=1:min(downtime_init[i],nPeriods)],
        UseF[i,t] == 0)
@addConstraint(m, FAST_MIN_DOWN_S1[i=FASTGENS, t=PERIODS, tt=max(1,t-downtime[i]+1):t],
        1 - UseF[i,t] >= DownF[i,tt])
@addConstraint(m, FAST_MIN_DOWN_S2[i=FASTGENS, t=PERIODS],
        1 - UseF[i,t] >= sum{DownF[i,tt], tt=max(1,t-downtime[i]+1):t})

# Min up time
@addConstraint(m, FAST_MIN_UP_INIT[i=FASTGENS, t=1:min(uptime_init[i],nPeriods)],
        UseF[i,t] == 1)
@addConstraint(m, FAST_MIN_UP_S1[i=FASTGENS, t=PERIODS, tt=max(1,t-uptime[i]+1):t],
        UseF[i,t] >= UpF[i,tt])
@addConstraint(m, FAST_MIN_UP_S2[i=FASTGENS, t=PERIODS],
        UseF[i,t] >= sum{UpF[i,tt], tt=max(1,t-uptime[i]+1):t})

# Ramping rate in normal operating status
@addConstraint(m, RAMP_DOWN0[i=GENERATORS],
        gen_0[i] - Gen[i,1] <= ramp_rate[i])
@addConstraint(m, RAMP_DOWN[i=GENERATORS, t=2:nPeriods],
        Gen[i,t-1] - Gen[i,t] <= ramp_rate[i])
@addConstraint(m, RAMP_UP0[i=GENERATORS],
        Gen[i,1] - gen_0[i] + Spin_Resv[i,1] <= ramp_rate[i])
@addConstraint(m, RAMP_UP[i=GENERATORS, t=2:nPeriods],
        Gen[i,t] - Gen[i,t-1] + Spin_Resv[i,t] <= ramp_rate[i]);

endCons7 = MathProgBase.numlinconstr(m);

# Spinning reserve requirement for system
@addConstraint(m, SPIN_RESV_REQ[t=PERIODS],
        sum{Spin_Resv[i,t], i=GENERATORS}
        >= spin_resv_rate * (total_demand[t] - total_wind_scen[t,1]));

# Spinning reserve capacity for individual unit
@addConstraint(m, SPIN_RESV_MAX_SLOW[i=SLOWGENS, t=PERIODS],
        Spin_Resv[i,t] <= spin_notice / 60. * ramp_rate[i] * Use[i,t])
@addConstraint(m, SPIN_RESV_MAX_FAST[i=FASTGENS, t=PERIODS],
        Spin_Resv[i,t] <= spin_notice / 60. * ramp_rate[i] * UseF[i,t])

# Power output capacity constraints
@addConstraint(m, POWER_OUTPUT_SLOW[i=SLOWGENS, t=PERIODS],
        Gen[i,t] == min_gen[i] * Use[i,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
@addConstraint(m, POWER_SEGMENT_SLOW[i=SLOWGENS, k=SEGMENTS, t=PERIODS],
        Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * Use[i,t])
@addConstraint(m, POWER_MAX_SLOW[i=SLOWGENS, t=PERIODS],
        Gen[i,t] + Spin_Resv[i,t] <= max_gen[i] * Use[i,t])
@addConstraint(m, POWER_OUTPUT_FAST[i=FASTGENS, t=PERIODS],
        Gen[i,t] == min_gen[i] * UseF[i,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
@addConstraint(m, POWER_SEGMENT_FAST[i=FASTGENS, k=SEGMENTS, t=PERIODS],
        Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * UseF[i,t])
@addConstraint(m, POWER_MAX_FAST[i=FASTGENS, t=PERIODS],
        Gen[i,t] + Spin_Resv[i,t] <= max_gen[i] * UseF[i,t])

endCons10 = MathProgBase.numlinconstr(m);

# Power balance constraints for system
@addConstraint(m, POWER_BALANCE[t=PERIODS],
        sum{Gen[i,t], i=GENERATORS} == total_demand[t] - total_wind_scen[t,1]);

endCons11 = MathProgBase.numlinconstr(m);

# Transmission constraints with load shift factor
@addConstraint(m, FLOW_BRANCH_LSF[l=BRANCHES, t=PERIODS],
	Flow[l,t] - sum{load_shift_factor[n,l] * Gen[i,t], n=BUSES, i=GENERATORS; gen_bus_id[i] == n}
        == sum{load_shift_factor[n,l] * wind_scen[wn,t,1], n=BUSES, wn=WINDS; wind_bus_id[wn] == n}
        - sum{load_shift_factor[n,l] * demand[n,t], n=BUSES});

# file name
fname = "../smps/ieee118-$nScenarios";

# write MPS file
println("Writing $fname.cor ...");
writeMPS(m, "$fname.cor");

# write .tim file
println("Writing $fname.tim ...");
fp = open("$fname.tim", "w");
print(fp,"*23456789 123456789 123456789 123456789 123456789 123456789\n");
print(fp,"TIME          JuMPModel\n");
print(fp,"PERIODS       IMPLICIT\n");
print(fp,"    VAR1      CON1                     PER1\n");
print(fp,"    VAR$startCol4      CON$startRow4                     PER2\n");
print(fp,"ENDATA");
close(fp);

# write .sto file
println("Writing $fname.sto ...");
prob = 1/nScenarios;
fp = open("$fname.sto", "w");
print(fp,"*23456789 123456789 123456789 123456789 123456789 123456789\n");
print(fp,"STOCH         JuMPModel\n");
print(fp,"SCENARIOS\n");
for s = 1:nScenarios
	print(fp, " SC SCEN$s     ROOT      $prob       PER2\n");
	for t = 1:length(PERIODS)
		i = endCons7 + t;
		rhs = spin_resv_rate * (total_demand[PERIODS[t]] - total_wind_scen[PERIODS[t],s]);
		print(fp, "    rhs    CON$i    $rhs\n");
	end
	for t = 1:length(PERIODS)
		i = endCons10 + t;
		rhs = total_demand[PERIODS[t]] - total_wind_scen[PERIODS[t],s];
		print(fp, "    rhs    CON$i    $rhs\n");
	end
	for l = 1:length(BRANCHES)
		for t = 1:length(PERIODS)
			i = endCons11 + (l-1)*length(PERIODS) + t;
			rhs = 0.0;
			for n = BUSES
				for wn = WINDS
					if wind_bus_id[wn] == n
						rhs += load_shift_factor[n,BRANCHES[l]] * wind_scen[wn,PERIODS[t],s];
					end
				end
				rhs -= load_shift_factor[n,BRANCHES[l]] * demand[n,PERIODS[t]];
			end
			print(fp, "    rhs    CON$i    $rhs\n");
		end
	end
end
print(fp,"ENDDATA");
close(fp);

