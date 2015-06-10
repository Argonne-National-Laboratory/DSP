
# Kibaek Kim - 2015 ANL MCS

# Check if JuMP is installed.
if Pkg.installed("JuMP") == nothing
	Pkg.add("JuMP")
end
using JuMP

const DSP_STAT_OPTIMAL          = 3000;
const DSP_STAT_PRIM_INFEASIBLE  = 3001;
const DSP_STAT_DUAL_INFEASIBLE  = 3002;
const DSP_STAT_LIM_ITERorTIME   = 3004;
const DSP_STAT_STOPPED_GAP      = 3005;
const DSP_STAT_STOPPED_NODE     = 3006;
const DSP_STAT_STOPPED_TIME     = 3007;
const DSP_STAT_STOPPED_USER     = 3008;
const DSP_STAT_STOPPED_SOLUTION = 3009;
const DSP_STAT_STOPPED_ITER     = 3010;
const DSP_STAT_STOPPED_UNKNOWN  = 3011;
const DSP_STAT_STOPPED_MPI      = 3012;
const DSP_STAT_ABORT            = 3013;
const DSP_STAT_LIM_PRIM_OBJ     = 3014;
const DSP_STAT_LIM_DUAL_OBJ     = 3015;
const DSP_STAT_UNKNOWN          = 3999;

macro dsp_ccall(func, args...)
	@unix_only return quote
		ccall(($func, "libDsp"), $(args...))
	end
	@windows_only return quote
		ccall(($func, "libDsp"), stdcall, $(args...))
	end
end

type DSP
	p::Ptr{Void}
	function DSP()
		p = @dsp_ccall("createEnv", Ptr{Void}, ())
		env = new(p)
		finalizer(env, freeDSP)
		return env
	end
end

type BranchingHyperplane
	nzcnt::Int
	indices::Vector{Int}
	values::Vector{Cdouble}
	priority::Int
end

function freeDSP(env::DSP)
	if env.p == C_NULL
		return
	end
	@dsp_ccall("freeEnv", Void, (Ptr{Void},), env.p)
	env.p = C_NULL
	return
end

function check_problem(env::DSP)
	if env.p == C_NULL
		error("Invalid DSP environment pointer")
	end
	return true
end

function freeModel(env::DSP)
	check_problem(env)
	@dsp_ccall("freeTssModel", Void, (Ptr{Void},), env.p)
end

function freeSolver(env::DSP)
	check_problem(env)
	@dsp_ccall("freeTssSolver", Void, (Ptr{Void},), env.p)
end

function getDataFormat(model::Model)
	# Get a column-wise sparse matrix
	mat = prepConstrMatrix(model)
	
	# Tranpose; now I have row-wise sparse matrix
	mat = mat'
	
	# sparse description
	start = convert(Vector{Cint}, mat.colptr - 1)
	index = convert(Vector{Cint}, mat.rowval - 1)
	value = mat.nzval
	
	# column type
	ctype = ""
	for i = 1:length(model.colCat)
		if model.colCat[i] == :Int
			ctype = ctype * "I";
		elseif model.colCat[i] == :Bin
			ctype = ctype * "B";
		else
			ctype = ctype * "C";
		end
	end
	ctype = convert(Vector{Uint8}, ctype)
	
	# objective coefficients
	obj, rlbd, rubd = JuMP.prepProblemBounds(model)
	
	return start, index, value, model.colLower, model.colUpper, ctype, obj, rlbd, rubd
end

function readSmps(env::DSP, filename)
	# Check pointer to TssModel
	check_problem(env)
	@dsp_ccall("readSmps", Void, (Ptr{Void}, Ptr{Uint8}), env.p, convert(Vector{Uint8}, filename))
end

function loadProblem(env::DSP, model::Model)
	# Check pointer to TssModel
	check_problem(env)
	# get scenario problem
	stoch = getStochastic(model)
	
	nscen  = convert(Cint, stoch.num_scen)
	ncols1 = convert(Cint, model.numCols)
	nrows1 = convert(Cint, length(model.linconstr))
	ncols2 = 0
	nrows2 = 0
	for s in 1:nscen
		if stoch.children[s].ext[:Skip] == false
			ncols2 = convert(Cint, stoch.children[s].numCols)
			nrows2 = convert(Cint, length(stoch.children[s].linconstr))
			break;
		end
	end
	
	#println("Number of scenarios: ", nscen)
	#println(" ncols1 ", ncols1, " nrows1 ", nrows1, " ncols2 ", ncols2, " nrows2 ", nrows2)
	@dsp_ccall("setNumberOfScenarios", Void, (Ptr{Void}, Cint), env.p, nscen)
	@dsp_ccall("setDimensions", Void, 
		(Ptr{Void}, Cint, Cint, Cint, Cint), 
		env.p, ncols1, nrows1, ncols2, nrows2)
	
	# get problem data
	start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(model)
	
	@dsp_ccall("loadFirstStage", Void, 
		(Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Uint8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
		env.p, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
	
	for s in 1:nscen
		# get model
		sb = stoch.children[s]
		if sb.ext[:Skip] == true
			#println("Skip ", s)
			continue;
		end
		probability = stoch.probability[s]
		# get model data
		start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(sb)
		@dsp_ccall("loadSecondStage", Void, 
			(Ptr{Void}, Cint, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Uint8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
			env.p, s-1, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
	end
	
end

function branchPriorities(env::DSP, priorities::Integer)
	@dsp_ccall("branchPriorities", Void, 
		(Ptr{Void}, Ptr{Cint}), env.p, convert(Vector{Cint}, priorities))
end

function addBranchingObjects(env::DSP, objects::Vector{BranchingHyperplane})
	for h in objects
		@dsp_ccall("addBranchingObject", Void, 
			(Ptr{Void}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint), 
			env.p, h.nzcnt, convert(Vector{Cint}, h.indices), convert(Vector{Cdouble}, h.values), h.priority)
	end
end

function evaluateSolution(env::DSP, solution)
	@dsp_ccall("evaluateSolution", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, convert(Vector{Cdouble}, solution))
end

function solveDe(env::DSP)
	@dsp_ccall("solveDe", Void, (Ptr{Void},), env.p)
end

function solveBd(env::DSP, nauxvars::Integer)
	@dsp_ccall("solveBd", Void, (Ptr{Void}, Cint), env.p, convert(Cint, nauxvars))
end

function setLogLevel(env::DSP, level::Integer)
	@dsp_ccall("setLogLevel", Void, (Ptr{Void}, Cint), env.p, convert(Cint, level))
end

function setNumCores(env::DSP, num::Integer)
	@dsp_ccall("setNumCores", Void, (Ptr{Void}, Cint), env.p, convert(Cint, num))
end

function setNodeLimit(env::DSP, num::Integer)
	@dsp_ccall("setNodeLimit", Void, (Ptr{Void}, Cint), env.p, convert(Cint, num))
end

function setIterLimit(env::DSP, num::Integer)
	@dsp_ccall("setIterLimit", Void, (Ptr{Void}, Cint), env.p, convert(Cint, num))
end

function setWallLimit(env::DSP, lim::Number)
	@dsp_ccall("setWallLimit", Void, (Ptr{Void}, Cdouble), env.p, convert(Cdouble, lim))
end

function setIntRelax(env::DSP, stage::Integer)
	@dsp_ccall("setIntRelax", Void, (Ptr{Void}, Cint), env.p, convert(Cint, stage))
end

function setBdAugScenarios(env::DSP, num::Integer, scenarios::Array{Int32,1})
	@dsp_ccall("setBdAugScenarios", Void, (Ptr{Void}, Cint, Ptr{Cint}), env.p, convert(Cint, num), scenarios)
end

function setBdAugScenarios(env::DSP, num::Integer, scenarios::Array{Int64,1})
	scenarios = convert(Vector{Cint}, scenarios)
	setBdAugScenarios(env, num, scenarios)
end

function setBendersAggressive(env::DSP, aggressive::Integer)
	@dsp_ccall("setBendersAggressive", Void, (Ptr{Void}, Cint), env.p, convert(Cint, aggressive))
end

function setDdAddFeasCuts(env::DSP, freq::Integer)
	@dsp_ccall("setDdAddFeasCuts", Void, (Ptr{Void}, Cint), env.p, convert(Cint, freq))
end

function setDdAddOptCuts(env::DSP, freq::Integer)
	@dsp_ccall("setDdAddOptCuts", Void, (Ptr{Void}, Cint), env.p, convert(Cint, freq))
end

function setDdEvalUb(env::DSP, freq::Integer)
	@dsp_ccall("setDdEvalUb", Void, (Ptr{Void}, Cint), env.p, convert(Cint, freq))
end

function setDdDualVarsLog(env::DSP, yesNo::Integer)
	@dsp_ccall("setDdDualVarsLog", Void, (Ptr{Void}, Cint), env.p, convert(Cint, yesNo))
end

function setDdCacheRecourse(env::DSP, yesNo::Integer)
	@dsp_ccall("setDdCacheRecourse", Void, (Ptr{Void}, Cint), env.p, convert(Cint, yesNo))
end

function setDdMasterSolver(env::DSP, solverType::Integer)
	@dsp_ccall("setDdMasterSolver", Void, (Ptr{Void}, Cint), env.p, convert(Cint, solverType))
end

function setDdMasterNumCutsPerIter(env::DSP, num::Integer)
	@dsp_ccall("setDdMasterNumCutsPerIter", Void, (Ptr{Void}, Cint), env.p, convert(Cint, num))
end

function setDdStoppingTolerance(env::DSP, tol::Number)
        @dsp_ccall("setDdStoppingTolerance", Void, (Ptr{Void}, Cdouble), env.p, convert(Cdouble, tol))
end

function setScipDisplayFreq(env::DSP, freq::Integer)
	@dsp_ccall("setScipDisplayFreq", Void, (Ptr{Void}, Cint), env.p, convert(Cint, freq))
end

function setScipLimitsGap(env::DSP, gap::Number)
	@dsp_ccall("setScipLimitsGap", Void, (Ptr{Void}, Cdouble), env.p, convert(Cdouble, gap))
end

function setScipLimitsTime(env::DSP, time::Number)
	@dsp_ccall("setScipLimitsTime", Void, (Ptr{Void}, Cdouble), env.p, convert(Cdouble, time))
end

function getNumScenarios(env::DSP)
	return @dsp_ccall("getNumScenarios", Cint, (Ptr{Void},), env.p)
end

function getNumRows(env::DSP, stage::Integer)
	return @dsp_ccall("getNumRows", Cint, (Ptr{Void}, Cint), env.p, convert(Cint, stage))
end

function getNumRows(env::DSP, stage::Integer)
	return @dsp_ccall("getNumRows", Cint, (Ptr{Void}, Cint), env.p, convert(Cint, stage))
end

function getNumCols(env::DSP, stage::Integer)
	return @dsp_ccall("getNumCols", Cint, (Ptr{Void}, Cint), env.p, convert(Cint, stage))
end

function getNumCols(env::DSP)
	num = getNumCols(env,0) + getNumScenarios(env) * getNumCols(env,1)
	return num
end

function getObjCoef(env::DSP)
	num = getNumCols(env)
	obj = Array(Cdouble, num)
	@dsp_ccall("getObjCoef", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, obj)
	return obj
end

function getSolutionTime(env::DSP)
	return @dsp_ccall("getSolutionTime", Cdouble, (Ptr{Void},), env.p)
end

function getSolutionStatus(env::DSP)
	status = @dsp_ccall("getSolutionStatus", Cint, (Ptr{Void},), env.p);
	if status == DSP_STAT_OPTIMAL
		return :Optimal
	elseif status == DSP_STAT_PRIM_INFEASIBLE
		return :PrimalInfeasible
	elseif status == DSP_STAT_DUAL_INFEASIBLE
		return :DualInfeasible
	elseif status == DSP_STAT_LIM_ITERorTIME
		return :IterOrTimeLimit
	elseif status == DSP_STAT_STOPPED_GAP
		return :StoppedGap
	elseif status == DSP_STAT_STOPPED_NODE
		return :StoppedNode
	elseif status == DSP_STAT_STOPPED_TIME
		return :StoppedTime
	elseif status == DSP_STAT_STOPPED_USER
		return :StoppedUser
	elseif status == DSP_STAT_STOPPED_SOLUTION
		return :StoppedSolution
	elseif status == DSP_STAT_STOPPED_ITER
		return :StoppedIter
	elseif status == DSP_STAT_STOPPED_UNKNOWN
		return :StoppedUnknown
	elseif status == DSP_STAT_STOPPED_MPI
		return :StoppedMPI
	elseif status == DSP_STAT_ABORT
		return :Abort
	elseif status == DSP_STAT_LIM_PRIM_OBJ
		return :PrimalObjLimit
	elseif status == DSP_STAT_LIM_DUAL_OBJ
		return :DualObjLimit
	else
		return :Unknown
	end
end

function getObjValue(env::DSP)
	return @dsp_ccall("getObjValue", Cdouble, (Ptr{Void},), env.p)
end

function getPrimalBound(env::DSP)
	return @dsp_ccall("getPrimalBound", Cdouble, (Ptr{Void},), env.p)
end

function getDualBound(env::DSP)
	return @dsp_ccall("getDualBound", Cdouble, (Ptr{Void},), env.p)
end

function getSolution(env::DSP, num::Integer)
	sol = Array(Cdouble, num)
	@dsp_ccall("getSolution", Void, (Ptr{Void}, Cint, Ptr{Cdouble}), env.p, num, sol)
	return sol
end

function getSolution(env::DSP)
	num = getNumCols(env)
	return getSolution(env,num)
end

function getNumIterations(env::DSP)
	return @dsp_ccall("getNumIterations", Cint, (Ptr{Void},), env.p)
end

function getNumNodes(env::DSP)
	return @dsp_ccall("getNumNodes", Cint, (Ptr{Void},), env.p)
end

function getDdNumInfeasSolutions(env::DSP)
	return @dsp_ccall("getDdNumInfeasSolutions", Cint, (Ptr{Void},), env.p)
end

function getDdIterTime(env::DSP)
	num = @dsp_ccall("getDdIterTimeSize", Cint, (Ptr{Void},), env.p)
	time = Array(Cdouble, num)
	@dsp_ccall("getDdIterTime", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, time)
	return time
end

function getDdMasterTime(env::DSP)
	num = @dsp_ccall("getDdMasterTimeSize", Cint, (Ptr{Void},), env.p)
	time = Array(Cdouble, num)
	@dsp_ccall("getDdMasterTime", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, time)
	return time
end

function getDdSubprobTime(env::DSP)
	num = @dsp_ccall("getDdSubprobTimeSize", Cint, (Ptr{Void},), env.p)
	time = Array(Cdouble, num)
	@dsp_ccall("getDdSubprobTime", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, time)
	return time
end

function getDdMasterObjValues(env::DSP)
	num = @dsp_ccall("getDdNumMasterObjValues", Cint, (Ptr{Void},), env.p)
	vals = Array(Cdouble, num)
	@dsp_ccall("getDdMasterObjValues", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, vals)
	return vals
end

function getDdSubproblemObjValues(env::DSP)
	num = @dsp_ccall("getDdNumSubproblemObjValues", Cint, (Ptr{Void},), env.p)
	vals = Array(Cdouble, num)
	@dsp_ccall("getDdSubproblemObjValues", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, vals)
	return vals
end

function getDdPrimalBounds(env::DSP)
	num = @dsp_ccall("getDdNumPrimalBounds", Cint, (Ptr{Void},), env.p)
	vals = Array(Cdouble, num)
	@dsp_ccall("getDdPrimalBounds", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, vals)
	return vals
end

function getDdDualBounds(env::DSP)
	num = @dsp_ccall("getDdNumDualBounds", Cint, (Ptr{Void},), env.p)
	vals = Array(Cdouble, num)
	@dsp_ccall("getDdDualBounds", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, vals)
	return vals
end

function getDdCpuTime(env::DSP)
	return @dsp_ccall("getDdCpuTime", Cdouble, (Ptr{Void},), env.p)
end

function getDdNumChangesOfMultiplier(env::DSP)
	return @dsp_ccall("getDdNumChangesOfMultiplier", Cint, (Ptr{Void},), env.p)
end

function getDdChangesOfMultiplier(env::DSP)
	num = getDdNumChangesOfMultiplier(env)
	changes = Array(Cdouble, num);
	@dsp_ccall("getDdChangesOfMultiplier", Void, (Ptr{Void}, Ptr{Cdouble}), env.p, changes)
	return changes
end

function printModel(env::DSP)
	@dsp_ccall("printModel", Void, (Ptr{Void},), env.p)
end

