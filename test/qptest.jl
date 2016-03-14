# Julia QP test

using JuMP, CPLEX

m = Model();

known = [-6.286111e+04
-3.032889e+04
-1.520000e+04
+4.011111e+01
+2.344444e+01
-6.355556e+01
-1.677778e+01
-2.811111e+01
+4.488889e+01
-2.333333e+01
+4.666667e+00
+1.866667e+01];

@defVar(m, x[i=0:11] == known[i+1]);

@setObjective(m, Max, sum{x[i], i=0:2});

#@addConstraint(m, CONS1[i=0:2], sum{x[j*3+i], j=1:3} == 0);

@addConstraint(m, CUT0,  x[0] - 183*x[3] - 67*x[4]   - 250*x[5]  <= -55883.3);
@addConstraint(m, CUT1,  x[1] - 120*x[6] - 80*x[7]   - 300*x[8]  <= -39533.3);
@addConstraint(m, CUT2,  x[2] - 100*x[9] - 25*x[10]  - 375*x[11] <= -19983.3);
@addConstraint(m, CUT3,  x[1] - 200*x[6]             - 300*x[8]  <= -36200);
@addConstraint(m, CUT4,  x[2] - 400*x[9] - 100*x[10]             <= 215708);
@addConstraint(m, CUT5,  x[0]                        - 500*x[5]  <= 177465);
if false
@addConstraint(m, CUT6,  x[1]                        - 500*x[8]  <= -184685);
@addConstraint(m, CUT7,  x[2]            - 500*x[10]             <= 25646.6);
@addConstraint(m, CUT8,  x[0]                        - 500*x[5]  <= -41588.2);
@addConstraint(m, CUT9,  x[1] - 420*x[6] - 80*x[7]               <= -21033.3);
@addConstraint(m, CUT10, x[1]            - 500*x[7]              <= 57403.4);
@addConstraint(m, CUT11, x[2]                        - 500*x[11] <= -26848.5);
@addConstraint(m, CUT12, x[0] - 500*x[3]                         <= -44792.4);
@addConstraint(m, CUT13, x[1]            - 500*x[7]              <= -11890.2);
@addConstraint(m, CUT14, x[0] - 250*x[3]             - 250*x[5]  <= -68125.5);
@addConstraint(m, CUT15, x[2]            - 125*x[10] - 375*x[11] <= -54148.9);
@addConstraint(m, CUT16, x[1]            - 500*x[7]              <= -19793.1);
@addConstraint(m, CUT17, x[0] - 250*x[3]             - 250*x[5]  <= -108733);
@addConstraint(m, CUT18, x[2]                        - 500*x[11] <= -29236.7);
end
solve(m)

#sum_of_thetas = sum(getValue(theta));
println(getObjectiveValue(m));
