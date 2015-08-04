# Julia script for unit commitment problem
# Kibaek Kim - 2015 ANL MCS

# ---------------
# Read data files
# ---------------

# limit_minupdown = false;
# nPeriods    = 12;
# nPeriods    = 24;
# nPeriods    = 48;
# nPeriods    = 96;

if isdefined(:nPeriods) == false
    println("Error: Number of periods not specified");
    exit();
end

if isdefined(:limit_minupdown) == false
    println("Error: Min up down limit not specified");
    exit();
end

tmpdat1 = readdlm("IEEE118/generator.dat", '\t');
tmpdat2 = readdlm("IEEE118/generator_cost_function.dat", '\t');
# tmpdat3 = readdlm("IEEE118/load_profile.dat", '\t');
tmpdat3 = readdlm("IEEE118/load_profile_$(nPeriods).dat", '\t');
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
# nPeriods    = 24;
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

wind_bus_id = tmpdat6[:,2]; # Wind power generator ID

# Load shift factor (This is pre-computed using branch_from_bus and branch_to_bus.)
load_shift_factor = tmpdat7';
flow_max          = tmpdat8[:,8]; # Transmission line capacity

# Min up/down limit
if limit_minupdown
    for i in 1:length(uptime)
            uptime[i] = min(uptime[i], 5)
    end
    for i in 1:length(downtime)
            downtime[i] = min(downtime[i], 5)
    end
end

# Scale problem according to number of periods
for i in 1:length(use_history)
        use_history[i] = ceil(use_history[i] * (nPeriods / 24.0))
end
for i in 1:length(uptime)
        uptime[i] = ceil(uptime[i] * (nPeriods / 24.0))
end
for i in 1:length(downtime)
        downtime[i] = ceil(downtime[i] * (nPeriods / 24.0))
end
# for i in 1:length(ramp_rate)
#         ramp_rate[i] = ceil(ramp_rate[i] / (nPeriods / 24.0))
# end
for i in 1:length(fixed_cost_gen)
        fixed_cost_gen[i] = ceil(fixed_cost_gen[i] / (nPeriods / 24.0))
end
for i in 1:length(cost_start)
        cost_start[i] = ceil(cost_start[i] / (nPeriods / 24.0))
end
for i in 1:length(cost_gen)
        cost_gen[i] = ceil(cost_gen[i] / (nPeriods / 24.0))
end

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
        downtime_init[i] = max(0, downtime[i] + use_history[i]) * (1 - use_0[i]);
        uptime_init[i]   = max(0, uptime[i] - use_history[i]) * use_0[i];
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
        winddat = readdlm(string("IEEE118/wp64-", nPeriods, "-", i, ".txt"));
        # winddat = readdlm(string("IEEE118/wp64-", i, ".txt"));
        wind_scen[i,:,:] = winddat[:,1:nScenarios];
end
total_wind_scen = reshape(sum(wind_scen,1), nPeriods, nScenarios);

# ----------------
# StochJuMP object
# ----------------
m = StochasticModel(nScenarios);

cc = DSPsolver.CouplingConstraints();

# ---------------------
# First-stage Variables
# ---------------------
# @defVar(m, Use[i=SLOWGENS, t=PERIODS], Bin)       # Generator on/off indicator
# @defVar(m, 0 <= Up[i=SLOWGENS, t=PERIODS] <= 1)   # Start up indicator
# @defVar(m, 0 <= Down[i=SLOWGENS, t=PERIODS] <= 1) # Shut down indicator

nIntervals = nPeriods / 6
assert(nPeriods % nIntervals == 0)
RANGES = {((i-1)*nPeriods/nIntervals+1):i*nPeriods/nIntervals for i in 1:nIntervals}
RANGES1 = vcat({2:nPeriods/nIntervals}, {((i-1)*nPeriods/nIntervals+1):i*nPeriods/nIntervals for i in 2:nIntervals})

INTERVALS = 1:nIntervals

# Main variables per interval
@defVar(m, Use[i=SLOWGENS, j=INTERVALS, t=(RANGES[j][1]-1):RANGES[j][end]], Bin)                            # Generator on/off indicator
# @defVar(m, 0 <= Use[i=SLOWGENS, j=INTERVALS, t=(RANGES[j][1]-1):RANGES[j][end]] <= 1)                            # Generator on/off indicator
@defVar(m, 0 <= Up[i=SLOWGENS, j=INTERVALS, t=max(1, RANGES[j][1]-uptime[i]+1):RANGES[j][end]] <= 1)        # Start up indicator
@defVar(m, 0 <= Down[i=SLOWGENS, j=INTERVALS, t=max(1, RANGES[j][1]-downtime[i]+1):RANGES[j][end]] <= 1)    # Shut down indicator

# Only required if using integer coupling constraints: set duplicate Up/Down variables to binary, as they are not implied to be binary
# for i in SLOWGENS
#     for j in INTERVALS
#         for t in max(1, RANGES[j][1]-uptime[i]+1):(RANGES[j][1]-1)
#             JuMP.setCategory(Up[i,j,t], :Bin)
#         end
#         for t in max(1, RANGES[j][1]-downtime[i]+1):(RANGES[j][1]-1)
#             JuMP.setCategory(Down[i,j,t], :Bin)
#         end
#     end
# end

# Equality constraints
for j in 2:nIntervals
    st = RANGES[j][1]  # initial step from interval

    # @addConstraint(m, USE_EQ[i=SLOWGENS, j],
    #     Use[i,j-1,st-1] == Use[i,j,st-1])
    # @addConstraint(m, UP_EQ[i=SLOWGENS, j, t=max(1,st-uptime[i]+1):(st-1)],
    #     Up[i,j-1,t] == Up[i,j,t])
    # @addConstraint(m, DOWN_EQ[i=SLOWGENS, j, t=max(1,st-downtime[i]+1):(st-1)],
    #     Down[i,j-1,t] == Down[i,j,t])

    # @addConstraint(m, UP_DOWN_EQ[i=SLOWGENS, j, t=max(1,st-uptime[i]+1):(st-1)],
    #     Up[i,j-1,t] + Down[i,j-1,t] == Up[i,j,t] + Down[i,j,t])
    # @addConstraint(m, UP_DOWN_EQ[i=SLOWGENS, j],
    #        sum{ (st - tt) * (Up[i,j-1,tt] + Down[i,j-1,tt]), tt=max(1,st-uptime[i]+1):(st-1)} 
    #     == sum{ (st - tt) * (Up[i,j,tt] + Down[i,j,tt]),     tt=max(1,st-uptime[i]+1):(st-1)})


    for i in SLOWGENS
        DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
            Use[i,j-1,st-1] == Use[i,j,st-1]));

        for t in max(1,st-uptime[i]+1):(st-1)
            DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
                Up[i,j-1,t] == Up[i,j,t]));
        end
        for t in max(1,st-downtime[i]+1):(st-1)
            DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
                Down[i,j-1,t] == Down[i,j,t]));
        end

        # TODO: Assuming uptime = downtime
        # for t in max(1,st-uptime[i]+1):(st-1)
        #     DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
        #         (Up[i,j-1,t] + Down[i,j-1,t]) == (Up[i,j,t] + Down[i,j,t])));
        # end

        # DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
        #        sum{ (st - tt) * (Up[i,j-1,tt] + Down[i,j-1,tt]), tt=max(1,st-uptime[i]+1):(st-1)} 
        #     == sum{ (st - tt) * (Up[i,j,tt] + Down[i,j,tt]),     tt=max(1,st-uptime[i]+1):(st-1)}));

    end
end

# ------------------------------
# First-stage Objective function
# ------------------------------
@setObjective(m, Min,
        sum{cost_start[i] * Up[i,j,t], i=SLOWGENS, j=INTERVALS, t=RANGES[j]}
        + sum{fixed_cost_gen[i] * Use[i,j,t], i=SLOWGENS, j=INTERVALS, t=RANGES[j]})


# -----------------------
# First-stage Constraints
# -----------------------

# # Linking Use / Up / Down variables
# @addConstraint(m, LINKING_SHUT_DOWN0[i=SLOWGENS],
#         Down[i,1] <= use_0[i])
# @addConstraint(m, LINKING_SHUT_DOWN[i=SLOWGENS, t=2:nPeriods],
#         Use[i,t-1] >= Down[i,t])
# @addConstraint(m, LINKING_START_UP0[i=SLOWGENS],
#         Up[i,1] <= 1 - use_0[i])
# @addConstraint(m, LINKING_START_UP[i=SLOWGENS, t=2:nPeriods],
#         1 - Use[i,t-1] >= Up[i,t])
# @addConstraint(m, LINKING_BOTH0[i=SLOWGENS],
#         Use[i,1] - use_0[i] == Up[i,1] - Down[i,1])
# @addConstraint(m, LINKING_BOTH[i=SLOWGENS, t=2:nPeriods],
#         Use[i,t] - Use[i,t-1] == Up[i,t] - Down[i,t])

# # Min down time
# @addConstraint(m, MIN_DOWN_INIT[i=SLOWGENS, t=1:min(downtime_init[i],nPeriods)],
#         Use[i,t] == 0)
# @addConstraint(m, MIN_DOWN_S1[i=SLOWGENS, t=PERIODS, s=max(1,t-downtime[i]+1):t],
#         1 - Use[i,t] >= Down[i,s])
# @addConstraint(m, MIN_DOWN_S2[i=SLOWGENS, t=PERIODS],
#         1 - Use[i,t] >= sum{Down[i,s], s=max(1,t-downtime[i]+1):t})

# # Min up time
# @addConstraint(m, MIN_UP_INIT[i=SLOWGENS, t=1:min(uptime_init[i],nPeriods)],
#         Use[i,t] == 1)
# @addConstraint(m, MIN_UP_S1[i=SLOWGENS, t=PERIODS, s=max(1,t-uptime[i]+1):t],
#         Use[i,t] >= Up[i,s])
# @addConstraint(m, MIN_UP_S2[i=SLOWGENS, t=PERIODS],
#         Use[i,t] >= sum{Up[i,s], s=max(1,t-uptime[i]+1):t})

# Linking Use / Up / Down variables
@addConstraint(m, LINKING_SHUT_DOWN0[i=SLOWGENS],
        Down[i,1,1] <= use_0[i])
@addConstraint(m, LINKING_SHUT_DOWN[i=SLOWGENS, j=INTERVALS, t=RANGES1[j]],
        Use[i,j,t-1] >= Down[i,j,t])
@addConstraint(m, LINKING_START_UP0[i=SLOWGENS],
        Up[i,1,1] <= 1 - use_0[i])
@addConstraint(m, LINKING_START_UP[i=SLOWGENS, j=INTERVALS, t=RANGES1[j]],
        1 - Use[i,j,t-1] >= Up[i,j,t])
@addConstraint(m, LINKING_BOTH0[i=SLOWGENS],
        Use[i,1,1] - use_0[i] == Up[i,1,1] - Down[i,1,1])
@addConstraint(m, LINKING_BOTH[i=SLOWGENS, j=INTERVALS, t=RANGES1[j]],
        Use[i,j,t] - Use[i,j,t-1] == Up[i,j,t] - Down[i,j,t])

# Min down time
@addConstraint(m, MIN_DOWN_INIT[i=SLOWGENS, j=filter(j -> downtime_init[i] >= RANGES[j][1], INTERVALS),
                t=1:min(downtime_init[i],RANGES[j][end])],
        Use[i,j,t] == 0)
@addConstraint(m, MIN_DOWN_S1[i=SLOWGENS, j=INTERVALS, t=RANGES[j], s=max(1,t-downtime[i]+1):t],
        1 - Use[i,j,t] >= Down[i,j,s])
@addConstraint(m, MIN_DOWN_S2[i=SLOWGENS, j=INTERVALS, t=RANGES[j]],
        1 - Use[i,j,t] >= sum{Down[i,j,s], s=max(1,t-downtime[i]+1):t})

# Min up time
@addConstraint(m, MIN_UP_INIT[i=SLOWGENS, j=filter(j -> uptime_init[i] >= RANGES[j][1], INTERVALS),
                t=1:min(uptime_init[i],RANGES[j][end])],
        Use[i,j,t] == 1)
@addConstraint(m, MIN_UP_S1[i=SLOWGENS, j=INTERVALS, t=RANGES[j], s=max(1,t-uptime[i]+1):t],
        Use[i,j,t] >= Up[i,j,s])
@addConstraint(m, MIN_UP_S2[i=SLOWGENS, j=INTERVALS, t=RANGES[j]],
        Use[i,j,t] >= sum{Up[i,j,s], s=max(1,t-uptime[i]+1):t})


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
        # @defVar(sb, UseF[i=FASTGENS, t=PERIODS], Bin)        # Generator on/off indicator
        # @defVar(sb, 0 <= UpF[i=FASTGENS, t=PERIODS] <= 1)    # Start up indicator
        # @defVar(sb, 0 <= DownF[i=FASTGENS, t=PERIODS] <= 1)  # Shut down indicator

        # Main variables per interval
        @defVar(sb, UseF[i=FASTGENS, j=INTERVALS, t=(RANGES[j][1]-1):RANGES[j][end]], Bin)                            # Generator on/off indicator
        # @defVar(sb, 0 <= UseF[i=FASTGENS, j=INTERVALS, t=(RANGES[j][1]-1):RANGES[j][end]] <= 1)                            # Generator on/off indicator
        @defVar(sb, 0 <= UpF[i=FASTGENS, j=INTERVALS, t=max(1, RANGES[j][1]-uptime[i]+1):RANGES[j][end]] <= 1)        # Start up indicator
        @defVar(sb, 0 <= DownF[i=FASTGENS, j=INTERVALS, t=max(1, RANGES[j][1]-downtime[i]+1):RANGES[j][end]] <= 1)    # Shut down indicator

        # Only required if using integer coupling constraints: set duplicate Up/Down variables to binary, as they are not implied to be binary
        # for i in FASTGENS
        #     for j in INTERVALS
        #         for t in max(1, RANGES[j][1]-uptime[i]+1):(RANGES[j][1]-1)
        #             JuMP.setCategory(UpF[i,j,t], :Bin)
        #         end
        #         for t in max(1, RANGES[j][1]-downtime[i]+1):(RANGES[j][1]-1)
        #             JuMP.setCategory(DownF[i,j,t], :Bin)
        #         end
        #     end
        # end

        @defVar(sb, 0 <= Gen[i=GENERATORS, t=PERIODS] <= max_gen[i])       # Power generation

        # @defVar(sb, 0 <= Gen[i=GENERATORS, j=INTERVALS, t=(RANGES[j][1]-1):RANGES[j][end]] <= max_gen[i])       # Power generation
        @defVar(sb, 0 <= Gen_Sgmt[i=GENERATORS, k=SEGMENTS, t=PERIODS] <= max_gen_sgmt[i,k])
        @defVar(sb, 0 <= Spin_Resv[i=GENERATORS, t=PERIODS] <= spin_notice / 60. * ramp_rate[i]) # Spinning reserve

        # Equality constraints
        for j in 2:nIntervals
            st = RANGES[j][1]  # initial step from interval

            # @addConstraint(sb, USE_EQ[i=FASTGENS, j],
            #     UseF[i,j-1,st-1] == UseF[i,j,st-1])
            # @addConstraint(sb, UP_EQ[i=FASTGENS, j, t=max(1,st-uptime[i]+1):(st-1)],
            #     UpF[i,j-1,t] == UpF[i,j,t])
            # @addConstraint(sb, DOWN_EQ[i=FASTGENS, j, t=max(1,st-downtime[i]+1):(st-1)],
            #     DownF[i,j-1,t] == DownF[i,j,t])
            # @addConstraint(sb, GEN_EQ[i=GENERATORS, j],
            #     Gen[i,j-1,st-1] == Gen[i,j,st-1])

            for i in FASTGENS
                DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
                    UseF[i,j-1,st-1] == UseF[i,j,st-1]));

                for t in max(1,st-uptime[i]+1):(st-1)
                    DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
                        UpF[i,j-1,t] == UpF[i,j,t]));
                end
                for t in max(1,st-downtime[i]+1):(st-1)
                    DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
                        DownF[i,j-1,t] == DownF[i,j,t]));
                end

                # TODO: Assuming uptime = downtime
                # for t in max(1,st-uptime[i]+1):(st-1)
                #     DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
                #         (UpF[i,j-1,t] + DownF[i,j-1,t]) == (UpF[i,j,t] + DownF[i,j,t])));
                # end

                # DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
                #        sum{ (st - tt) * (UpF[i,j-1,tt] + DownF[i,j-1,tt]), tt=max(1,st-uptime[i]+1):(st-1)} 
                #     == sum{ (st - tt) * (UpF[i,j,tt] + DownF[i,j,tt]),     tt=max(1,st-uptime[i]+1):(st-1)}));
            end

            # for i in GENERATORS
            #     DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(
            #         Gen[i,j-1,st-1] == Gen[i,j,st-1]));
            # end
        end

        # -------------------------------
        # Second-stage Objective function
        # -------------------------------

        @setObjective(sb, Min,
                sum{cost_start[i] * UpF[i,j,t], i=FASTGENS, j=INTERVALS, t=RANGES[j]}
                + sum{fixed_cost_gen[i] * UseF[i,j,t], i=FASTGENS, j=INTERVALS, t=RANGES[j]}
                + sum{cost_gen[i,k] * Gen_Sgmt[i,k,t], i=GENERATORS, k=SEGMENTS, t=PERIODS})

        # ------------------------
        # Second-stage Constraints
        # ------------------------

        # # Linking Use / Up / Down variables
        # @addConstraint(sb, FAST_LINKING_SHUT_DOWN0[i=FASTGENS],
        #         DownF[i,1] <= use_0[i])
        # @addConstraint(sb, FAST_LINKING_SHUT_DOWN[i=FASTGENS, t=2:nPeriods],
        #         UseF[i,t-1] >= DownF[i,t])
        # @addConstraint(sb, FAST_LINKING_START_UP0[i=FASTGENS],
        #         UpF[i,1] <= 1 - use_0[i])
        # @addConstraint(sb, FAST_LINKING_START_UP[i=FASTGENS, t=2:nPeriods],
        #         1 - UseF[i,t-1] >= UpF[i,t])
        # @addConstraint(sb, FAST_LINKING_BOTH0[i=FASTGENS],
        #         UseF[i,1] - use_0[i] == UpF[i,1] - DownF[i,1])
        # @addConstraint(sb, FAST_LINKING_BOTH[i=FASTGENS, t=2:nPeriods],
        #         UseF[i,t] - UseF[i,t-1] == UpF[i,t] - DownF[i,t])

        # # Min down time
        # @addConstraint(sb, FAST_MIN_DOWN_INIT[i=FASTGENS, t=1:min(downtime_init[i],nPeriods)],
        #         UseF[i,t] == 0)
        # @addConstraint(sb, FAST_MIN_DOWN_S1[i=FASTGENS, t=PERIODS, tt=max(1,t-downtime[i]+1):t],
        #         1 - UseF[i,t] >= DownF[i,tt])
        # @addConstraint(sb, FAST_MIN_DOWN_S2[i=FASTGENS, t=PERIODS],
        #         1 - UseF[i,t] >= sum{DownF[i,tt], tt=max(1,t-downtime[i]+1):t})

        # # Min up time
        # @addConstraint(sb, FAST_MIN_UP_INIT[i=FASTGENS, t=1:min(uptime_init[i],nPeriods)],
        #         UseF[i,t] == 1)
        # @addConstraint(sb, FAST_MIN_UP_S1[i=FASTGENS, t=PERIODS, tt=max(1,t-uptime[i]+1):t],
        #         UseF[i,t] >= UpF[i,tt])
        # @addConstraint(sb, FAST_MIN_UP_S2[i=FASTGENS, t=PERIODS],
        #         UseF[i,t] >= sum{UpF[i,tt], tt=max(1,t-uptime[i]+1):t})

        # Linking Use / Up / Down variables
        @addConstraint(sb, FAST_LINKING_SHUT_DOWN0[i=FASTGENS],
                DownF[i,1,1] <= use_0[i])
        @addConstraint(sb, FAST_LINKING_SHUT_DOWN[i=FASTGENS, j=INTERVALS, t=RANGES1[j]],
                UseF[i,j,t-1] >= DownF[i,j,t])
        @addConstraint(sb, FAST_LINKING_START_UP0[i=FASTGENS],
                UpF[i,1,1] <= 1 - use_0[i])
        @addConstraint(sb, FAST_LINKING_START_UP[i=FASTGENS, j=INTERVALS, t=RANGES1[j]],
                1 - UseF[i,j,t-1] >= UpF[i,j,t])
        @addConstraint(sb, FAST_LINKING_BOTH0[i=FASTGENS],
                UseF[i,1,1] - use_0[i] == UpF[i,1,1] - DownF[i,1,1])
        @addConstraint(sb, FAST_LINKING_BOTH[i=FASTGENS, j=INTERVALS, t=RANGES1[j]],
                UseF[i,j,t] - UseF[i,j,t-1] == UpF[i,j,t] - DownF[i,j,t])

        # Min down time
        @addConstraint(sb, FAST_MIN_DOWN_INIT[i=FASTGENS, j=filter(j -> downtime_init[i] >= RANGES[j][1], INTERVALS),
                        t=1:min(downtime_init[i],RANGES[j][end])],
                UseF[i,j,t] == 0)
        @addConstraint(sb, FAST_MIN_DOWN_S1[i=FASTGENS, j=INTERVALS, t=RANGES[j], tt=max(1,t-downtime[i]+1):t],
                1 - UseF[i,j,t] >= DownF[i,j,tt])
        @addConstraint(sb, FAST_MIN_DOWN_S2[i=FASTGENS, j=INTERVALS, t=RANGES[j]],
                1 - UseF[i,j,t] >= sum{DownF[i,j,tt], tt=max(1,t-downtime[i]+1):t})

        # Min up time
        @addConstraint(sb, FAST_MIN_UP_INIT[i=FASTGENS, j=filter(j -> uptime_init[i] >= RANGES[j][1], INTERVALS),
                        t=1:min(uptime_init[i],RANGES[j][end])],
                UseF[i,j,t] == 1)
        @addConstraint(sb, FAST_MIN_UP_S1[i=FASTGENS, j=INTERVALS, t=RANGES[j], tt=max(1,t-uptime[i]+1):t],
                UseF[i,j,t] >= UpF[i,j,tt])
        @addConstraint(sb, FAST_MIN_UP_S2[i=FASTGENS, j=INTERVALS, t=RANGES[j]],
                UseF[i,j,t] >= sum{UpF[i,j,tt], tt=max(1,t-uptime[i]+1):t})


        # # Ramping rate in normal operating status
        # @addConstraint(sb, RAMP_DOWN0[i=GENERATORS],
        #         gen_0[i] - Gen[i,1,1] <= ramp_rate[i])
        # @addConstraint(sb, RAMP_DOWN[i=GENERATORS, j=INTERVALS, t=RANGES1[j]],
        #         Gen[i,j,t-1] - Gen[i,j,t] <= ramp_rate[i])
        # @addConstraint(sb, RAMP_UP0[i=GENERATORS],
        #         Gen[i,1,1] - gen_0[i] + Spin_Resv[i,1] <= ramp_rate[i])
        # @addConstraint(sb, RAMP_UP[i=GENERATORS, j=INTERVALS, t=RANGES1[j]],
        #         Gen[i,j,t] - Gen[i,j,t-1] + Spin_Resv[i,t] <= ramp_rate[i])

        # Ramping rate in normal operating status
        @addConstraint(sb, RAMP_DOWN0[i=GENERATORS],
                gen_0[i] - Gen[i,1] <= ramp_rate[i])
        for i in GENERATORS
            for t in 2:nPeriods
                DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(Gen[i,t-1] - Gen[i,t] <= ramp_rate[i]))
            end
        end
        @addConstraint(sb, RAMP_UP0[i=GENERATORS],
                Gen[i,1] - gen_0[i] + Spin_Resv[i,1] <= ramp_rate[i])
        for i in GENERATORS
            for t in 2:nPeriods
                DSPsolver.addCouplingConstraint(cc, @JuMP.LinearConstraint(Gen[i,t] - Gen[i,t-1] + Spin_Resv[i,t] <= ramp_rate[i]))
            end
        end

        # Spinning reserve requirement for system
        @addConstraint(sb, SPIN_RESV_REQ[t=PERIODS],
                sum{Spin_Resv[i,t], i=GENERATORS}
                >= spin_resv_rate * (total_demand[t] - total_wind_scen[t,s]))

        # Spinning reserve capacity for individual unit
        @addConstraint(sb, SPIN_RESV_MAX_SLOW[i=SLOWGENS, j=INTERVALS, t=RANGES[j]],
                Spin_Resv[i,t] <= spin_notice / 60. * ramp_rate[i] * Use[i,j,t])
        @addConstraint(sb, SPIN_RESV_MAX_FAST[i=FASTGENS, j=INTERVALS, t=RANGES[j]],
                Spin_Resv[i,t] <= spin_notice / 60. * ramp_rate[i] * UseF[i,j,t])

        # # Power output capacity constraints
        # @addConstraint(sb, POWER_OUTPUT_SLOW[i=SLOWGENS, j=INTERVALS, t=RANGES[j]],
        #         Gen[i,j,t] == min_gen[i] * Use[i,j,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
        # @addConstraint(sb, POWER_SEGMENT_SLOW[i=SLOWGENS, k=SEGMENTS, j=INTERVALS, t=RANGES[j]],
        #         Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * Use[i,j,t])
        # @addConstraint(sb, POWER_MAX_SLOW[i=SLOWGENS, j=INTERVALS, t=RANGES[j]],
        #         Gen[i,j,t] + Spin_Resv[i,t] <= max_gen[i] * Use[i,j,t])
        # @addConstraint(sb, POWER_OUTPUT_FAST[i=FASTGENS, j=INTERVALS, t=RANGES[j]],
        #         Gen[i,j,t] == min_gen[i] * UseF[i,j,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
        # @addConstraint(sb, POWER_SEGMENT_FAST[i=FASTGENS, k=SEGMENTS, j=INTERVALS, t=RANGES[j]],
        #         Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * UseF[i,j,t])
        # @addConstraint(sb, POWER_MAX_FAST[i=FASTGENS, j=INTERVALS, t=RANGES[j]],
        #         Gen[i,j,t] + Spin_Resv[i,t] <= max_gen[i] * UseF[i,j,t])

        # # Power balance constraints for system
        # @addConstraint(sb, POWER_BALANCE[j=INTERVALS, t=RANGES[j]],
        #         sum{Gen[i,j,t], i=GENERATORS} == total_demand[t] - total_wind_scen[t,s])

        # # Transmission constraints with load shift factor (These can be lazy constraints.)
        # @addConstraint(sb, FLOW_BRANCH_LSF_LB[l=BRANCHES, j=INTERVALS, t=RANGES[j]],
        #         sum{load_shift_factor[n,l] * Gen[i,j,t], n=BUSES, i=GENERATORS; gen_bus_id[i] == n}
        #         >= sum{load_shift_factor[n,l] * demand[n,t], n=BUSES}
        #         - sum{load_shift_factor[n,l] * wind_scen[wn,t,s], n=BUSES, wn=WINDS; wind_bus_id[wn] == n}
        #         - flow_max[l])
        # @addConstraint(sb, FLOW_BRANCH_LSF_UB[l=BRANCHES, j=INTERVALS, t=RANGES[j]],
        #         sum{load_shift_factor[n,l] * Gen[i,j,t], n=BUSES, i=GENERATORS; gen_bus_id[i] == n}
        #         <= sum{load_shift_factor[n,l] * demand[n,t], n=BUSES}
        #         - sum{load_shift_factor[n,l] * wind_scen[wn,t,s], n=BUSES, wn=WINDS; wind_bus_id[wn] == n}
        #         + flow_max[l])

        @addConstraint(sb, POWER_OUTPUT_SLOW[i=SLOWGENS, j=INTERVALS, t=RANGES[j]],
                Gen[i,t] == min_gen[i] * Use[i,j,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
        @addConstraint(sb, POWER_SEGMENT_SLOW[i=SLOWGENS, k=SEGMENTS, j=INTERVALS, t=RANGES[j]],
                Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * Use[i,j,t])
        @addConstraint(sb, POWER_MAX_SLOW[i=SLOWGENS, j=INTERVALS, t=RANGES[j]],
                Gen[i,t] + Spin_Resv[i,t] <= max_gen[i] * Use[i,j,t])
        @addConstraint(sb, POWER_OUTPUT_FAST[i=FASTGENS, j=INTERVALS, t=RANGES[j]],
                Gen[i,t] == min_gen[i] * UseF[i,j,t] + sum{Gen_Sgmt[i,k,t], k=SEGMENTS})
        @addConstraint(sb, POWER_SEGMENT_FAST[i=FASTGENS, k=SEGMENTS, j=INTERVALS, t=RANGES[j]],
                Gen_Sgmt[i,k,t] <= max_gen_sgmt[i,k] * UseF[i,j,t])
        @addConstraint(sb, POWER_MAX_FAST[i=FASTGENS, j=INTERVALS, t=RANGES[j]],
                Gen[i,t] + Spin_Resv[i,t] <= max_gen[i] * UseF[i,j,t])

        # Power balance constraints for system
        @addConstraint(sb, POWER_BALANCE[t=PERIODS],
                sum{Gen[i,t], i=GENERATORS} == total_demand[t] - total_wind_scen[t,s])

        # Transmission constraints with load shift factor (These can be lazy constraints.)
        @addConstraint(sb, FLOW_BRANCH_LSF_LB[l=BRANCHES, t=PERIODS],
                sum{load_shift_factor[n,l] * Gen[i,t], n=BUSES, i=GENERATORS; gen_bus_id[i] == n}
                >= sum{load_shift_factor[n,l] * demand[n,t], n=BUSES}
                - sum{load_shift_factor[n,l] * wind_scen[wn,t,s], n=BUSES, wn=WINDS; wind_bus_id[wn] == n}
                - flow_max[l])
        @addConstraint(sb, FLOW_BRANCH_LSF_UB[l=BRANCHES, t=PERIODS],
                sum{load_shift_factor[n,l] * Gen[i,t], n=BUSES, i=GENERATORS; gen_bus_id[i] == n}
                <= sum{load_shift_factor[n,l] * demand[n,t], n=BUSES}
                - sum{load_shift_factor[n,l] * wind_scen[wn,t,s], n=BUSES, wn=WINDS; wind_bus_id[wn] == n}
                + flow_max[l])
end

DSPsolver.loadProblem(m);

varBlocks = (Variable => Int)[];

for i in SLOWGENS, j in INTERVALS, t in (RANGES[j][1]-1):RANGES[j][end]
        varBlocks[Use[i,j,t]] = j;
end
for i in SLOWGENS, j in INTERVALS, t in max(1, RANGES[j][1]-uptime[i]+1):RANGES[j][end]
        varBlocks[Up[i,j,t]] = j;
end
for i in SLOWGENS, j in INTERVALS, t in max(1, RANGES[j][1]-downtime[i]+1):RANGES[j][end]
        varBlocks[Down[i,j,t]] = j;
end
for sb in getchildren(m)
        UseF = JuMP.getVar(sb, :UseF);
        UpF = JuMP.getVar(sb, :UpF);
        DownF = JuMP.getVar(sb, :DownF);
        Gen = JuMP.getVar(sb, :Gen);
        Gen_Sgmt = JuMP.getVar(sb, :Gen_Sgmt);
        Spin_Resv = JuMP.getVar(sb, :Spin_Resv);

        for i in FASTGENS, j in INTERVALS, t in (RANGES[j][1]-1):RANGES[j][end]
                varBlocks[UseF[i,j,t]] = j
        end
        for i in FASTGENS, j in INTERVALS, t in max(1, RANGES[j][1]-uptime[i]+1):RANGES[j][end]
                varBlocks[UpF[i,j,t]] = j
        end
        for i in FASTGENS, j in INTERVALS, t in max(1, RANGES[j][1]-downtime[i]+1):RANGES[j][end]
                varBlocks[DownF[i,j,t]] = j
        end

        # for i in GENERATORS, j in INTERVALS, t in (RANGES[j][1]-1):RANGES[j][end]
        #         varBlocks[Gen[i,j,t]] = j;
        # end
        for i in GENERATORS, j in INTERVALS, t in RANGES[j]
                varBlocks[Gen[i,t]] = j;
        end
        for i in GENERATORS, j in INTERVALS, t in RANGES[j]
                varBlocks[Spin_Resv[i,t]] = j
        end
        for i in GENERATORS, k in SEGMENTS, j in INTERVALS, t in RANGES[j]
                varBlocks[Gen_Sgmt[i,k,t]] = j
        end
end

DSPsolver.addVarBlocks(cc, varBlocks);

DSPsolver.loadDecomposition(m, cc);
