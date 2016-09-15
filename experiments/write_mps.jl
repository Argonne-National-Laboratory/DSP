# Julia script of writing MPS files from SMPS files
# Kibaek Kim - ANL MCS 2016

problems=[
	"dcap233_200"
	"dcap332_200"
	"dcap332_300"
	"dcap342_200"
	"dcap342_300"
	"sslp_5_25_50"
	"sslp_5_25_100"
	"sslp_5_50_50"
	"sslp_5_50_100"
	"sslp_10_50_50"
	"sslp_10_50_100"
	"sslp_15_45_5"
	"dcap233_300"
	"dcap233_500"
	"dcap243_200"
	"dcap243_300"
	"dcap243_500"
	"dcap332_500"
	"dcap342_500"
	"sslp_10_50_500"
	"sslp_10_50_1000"
	"sslp_10_50_2000"
	"sslp_15_45_10"
	"sslp_15_45_15"
]

using DSPsolver

for p in problems
	DSPsolver.readSmps("../smps/$p");
	DSPsolver.writeMps("./mps/$p.mps");
end

