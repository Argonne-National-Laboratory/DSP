# Kibaek Kim - ANL MCS 2014
# Farmer example from Birge and Louveaux book.


# STOCHASTIC MODELING FRAMEWORK

NS = 3;                        # number of scenarios
probability = [1/3, 1/3, 1/3]; # probability


# FIRST-STAGE MODEL

CROPS = 1:3; # set of crops (wheat, corn and sugar beets, resp.)
Cost = [150 230 260]; # cost of planting crops
Budget = 500; # budget capacity


# SECOND-STAGE MODELS

PURCH = 1:2; # set of crops to purchase (wheat and corn, resp.)
SELL  = 1:4; # set of crops to sell (wheat, corn, sugar beets under 6K and those over 6K)
Purchase = [238 210;
            238 210;
            238 210];   # purchase price
Sell = [170 150 36 10;
        170 150 36 10;
        170 150 36 10]; # selling price
Yield = [3.0 3.6 24.0;
         2.5 3.0 20.0;
         2.0 2.4 16.0];
Minreq = [200 240 0;
          200 240 0;
          200 240 0]; # minimum crop requirement

