# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.


# CREATE STOCHASTIC MODEL

m = StochasticModel(NS);


# FIRST-STAGE MODEL

# first-stage variables
@defVar(m, x[i=CROPS] >= 0, Int)

# first-stage objective
@setObjective(m, Min, sum{Cost[i] * x[i], i=CROPS})

# first-stage constraint
@addConstraint(m, const_budget,
               sum{x[i], i=CROPS} <= Budget)


# SECOND-STAGE MODELS

@second_stage m s begin
    # stochastic block
    sb = StochasticBlock(m, probability[s]);

    # second-stage variables
    @defVar(sb, y[j=PURCH] >= 0)
    @defVar(sb, w[k=SELL] >= 0)

    # objective
    @setObjective(sb, Min,
                  sum{Purchase[s,j] * y[j], j=PURCH}
                  - sum{Sell[s,k] * w[k], k=SELL})

    # constraints
    @addConstraint(sb, const_minreq[j=PURCH],
                   Yield[s,j] * x[j] + y[j] - w[j] >= Minreq[s,j])
    @addConstraint(sb, const_minreq_beets,
                   Yield[s,3] * x[3] - w[3] - w[4] >= Minreq[s,3])
    @addConstraint(sb, const_aux, w[3] <= 6000)
end


