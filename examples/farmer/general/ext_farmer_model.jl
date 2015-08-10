# ctjandra - ANL MCS 2015
# Farmer example from Birge and Louveaux book.

# This is an extensive form of the model from examples/julia/farmer with nonanticipativity
# constraints, written in a general decomposition format **for illustration purposes only**.
# DSP contains features specific to stochastic programs that will not be applied if the
# model is written in this form (as opposed to the form in examples/julia/farmer).

SCENARIOS = 1:NS

# CREATE MODEL

m = Model();

# MODEL IN EXTENSIVE FORM WITH NONANTICIPATIVITY CONSTRAINTS

# first-stage variables with copies for each stage
@defVar(m, x[s=SCENARIOS, i=CROPS] >= 0, Int)

# second-stage variables
@defVar(m, y[s=SCENARIOS, j=PURCH] >= 0)
@defVar(m, w[s=SCENARIOS, k=SELL] >= 0)

# first-stage constraint
@addConstraint(m, const_budget[s=SCENARIOS],
               sum{x[s,i], i=CROPS} <= Budget)

# second-stage constraints
@addConstraint(m, const_minreq[s=SCENARIOS, j=PURCH],
               Yield[s,j] * x[s,j] + y[s,j] - w[s,j] >= Minreq[s,j])
@addConstraint(m, const_minreq_beets[s=SCENARIOS],
               Yield[s,3] * x[s,3] - w[s,3] - w[s,4] >= Minreq[s,3])
@addConstraint(m, const_aux[s=SCENARIOS], w[s,3] <= 6000)

# objective
@setObjective(m, Min,
  sum{probability[s]
    * (sum{Cost[i] * x[s,i], i=CROPS}
      + sum{Purchase[s,j] * y[s,j], j=PURCH}
      - sum{Sell[s,k] * w[s,k], k=SELL}),
    s=SCENARIOS})

# nonanticipativity constraints (coupling constraints)
for s in 1:NS-1, i in CROPS
  DSPsolver.addCouplingConstraint(m, @LinearConstraint(x[s,i] == x[s+1,i]))
end

# variable partition
for s in SCENARIOS, i in CROPS
  DSPsolver.setVarSubproblem(m, x[s,i], s)
end
for s in SCENARIOS, j in PURCH
  DSPsolver.setVarSubproblem(m, y[s,j], s)
end
for s in SCENARIOS, k in SELL
  DSPsolver.setVarSubproblem(m, w[s,k], s)
end

