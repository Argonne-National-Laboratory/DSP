# Kibaek Kim - ANL MCS 2016
# Farmer example from Birge and Louveaux book.

using JuMP, Dsp

NS = 3;                        # number of scenarios
probability = [1/3, 1/3, 1/3]; # probability

CROPS  = 1:3 # set of crops (wheat, corn and sugar beets, resp.)
PURCH  = 1:2 # set of crops to purchase (wheat and corn, resp.)
SELL   = 1:4 # set of crops to sell (wheat, corn, sugar beets under 6K and those over 6K)

Cost     = [150 230 260]   # cost of planting crops
Budget   = 500             # budget capacity
Purchase = [238 210];      # purchase price
Sell     = [170 150 36 10] # selling price
Yield    = [3.0 3.6 24.0;
            2.5 3.0 20.0;
            2.0 2.4 16.0]
Minreq   = [200 240 0]     # minimum crop requirement

# JuMP model
m = Model(NS)

@variable(m, x[i=CROPS] >= 0, Int)
@objective(m, Min, sum(Cost[i] * x[i] for i=CROPS))
@constraint(m, const_budget, sum(x[i] for i=CROPS) <= Budget)

for s in blockids()
    blk = Model(m, s, probability[s])

    @variable(blk, y[j=PURCH] >= 0)
    @variable(blk, w[k=SELL] >= 0)

    @objective(blk, Min, sum(Purchase[j] * y[j] for j=PURCH) - sum(Sell[k] * w[k] for k=SELL))

    @constraint(blk, const_minreq[j=PURCH], Yield[s,j] * x[j] + y[j] - w[j] >= Minreq[j])
    @constraint(blk, const_minreq_beets, Yield[s,3] * x[3] - w[3] - w[4] >= Minreq[3])
    @constraint(blk, const_aux, w[3] <= 6000)
end
