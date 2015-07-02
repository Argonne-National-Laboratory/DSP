/*
 * SolverInterfaceScip.cpp
 *
 *  Created on: Dec 8, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** DSP */
#include "SolverInterface/SolverInterfaceScip.h"
#include "SolverInterface/SCIPconshdlrBenders.h"
#include "SolverInterface/SCIPconshdlrBendersDd.h"
#include "Utility/StoMacros.h"
#include "Utility/StoMessage.h"

/** SCIP */
#include "scip/scipdefplugins.h"
#include "../examples/Queens/src/scip_exception.hpp"

SolverInterfaceScip::SolverInterfaceScip(StoParam * par) :
	SolverInterface(par),
	scip_(NULL),
	solution_(NULL),
	nvars_(0),
	vars_(NULL),
	clbd_(NULL),
	cubd_(NULL),
	nconss_(0),
	conss_(NULL),
	rlbd_(NULL),
	rubd_(NULL)
{
	STO_RTN_CHECK_THROW(initialize(), "initialize", "SolverInterfaceScip");
}

/** copy constructor */
SolverInterfaceScip::SolverInterfaceScip(SolverInterfaceScip * si) :
	SolverInterface(si->par_),
	scip_(NULL)
{
	STO_RTN_CHECK_THROW(initialize(), "initialize", "SolverInterfaceScip");

	SCIP * source = si->getSCIP();

	/** copy problem */
	SCIP_CALL_ABORT(SCIPcopyOrigProb(source, scip_, NULL, NULL, "clone"));

	/** retrieve the number of variables */
	int nrows = SCIPgetNOrigConss(scip_);
	int ncols = SCIPgetNOrigVars(scip_);

	/** copy solution */
	if (si->getSolution())
	{
		solution_ = new double [ncols];
		CoinCopyN(si->getSolution(), ncols, solution_);
	}

	/** copy variables */
	SCIP_VAR ** vars = SCIPgetOrigVars(scip_);
	vars_ = new SCIP_VAR * [ncols];
	clbd_ = new double [ncols];
	cubd_ = new double [ncols];
	for (int j = 0; j < ncols; ++j)
	{
		SCIP_CALL_ABORT(SCIPcaptureVar(scip_, vars[j]));
		vars_[j] = vars[j];
	}
	CoinCopyN(si->getColLower(), ncols, clbd_);
	CoinCopyN(si->getColUpper(), ncols, cubd_);

	/** copy constraints */
	SCIP_CONS ** conss = SCIPgetOrigConss(scip_);
	conss_ = new SCIP_CONS * [nrows];
	rlbd_ = new double [nrows];
	rubd_ = new double [nrows];
	for (int i = 0; i < nrows; ++i)
	{
		SCIP_CALL_ABORT(SCIPcaptureCons(scip_, conss[i]));
		conss_[i] = conss[i];
	}
	CoinCopyN(si->getRowLower(), ncols, rlbd_);
	CoinCopyN(si->getRowUpper(), ncols, rubd_);
}

/** clone */
SolverInterface * SolverInterfaceScip::clone()
{
	return new SolverInterfaceScip(this);
}

SolverInterfaceScip::~SolverInterfaceScip()
{
	STO_RTN_CHECK_THROW(finalize(), "finalize", "SolverInterfaceScip");
	FREE_ARRAY_PTR(solution_);
	FREE_ARRAY_PTR(vars_);
	FREE_ARRAY_PTR(clbd_);
	FREE_ARRAY_PTR(cubd_);
	FREE_ARRAY_PTR(conss_);
	FREE_ARRAY_PTR(rlbd_);
	FREE_ARRAY_PTR(rubd_);
}

/** create SCIP */
STO_RTN_CODE SolverInterfaceScip::initialize()
{
	if (scip_ == NULL)
	{
		SCIP_CALL_ABORT(SCIPcreate(&scip_));
		SCIP_CALL_ABORT(SCIPincludeDefaultPlugins(scip_));
		setClockType(2);
		SCIP_CALL_ABORT(SCIPsetIntParam(scip_, "display/freq", par_->ScipDisplayFreq_));
		SCIP_CALL_ABORT(SCIPsetRealParam(scip_, "limits/gap", par_->ScipLimitsGap_));
	}
	return STO_RTN_OK;
}

/** free SCIP */
STO_RTN_CODE SolverInterfaceScip::finalize()
{
	if (scip_)
	{
		for (int j = 0; j < nvars_; ++j)
			SCIP_CALL_ABORT(SCIPreleaseVar(scip_, &vars_[j]));
		for (int i = 0; i < nconss_; ++i)
			SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &conss_[i]));
		SCIP_CALL_ABORT(SCIPfree(&scip_));
		nvars_ = 0;
		nconss_ = 0;
		scip_ = NULL;
	}
	return STO_RTN_OK;
}

/** load problem */
void SolverInterfaceScip::loadProblem(
		OsiSolverInterface * si,
		const char * probname)
{
	CoinPackedMatrix * mat = new CoinPackedMatrix(*si->getMatrixByRow());
	char * ctype = new char [si->getNumCols()];
	for (int j = 0; j < si->getNumCols(); ++j)
	{
		if (si->getColType()[j] == 0)
			ctype[j] = 'C';
		else if (si->getColType()[j] == 1)
			ctype[j] = 'B';
		else
			ctype[j] = 'I';
	}
	loadProblem(mat, si->getColLower(), si->getColUpper(), si->getObjCoefficients(),
			const_cast<const char*>(ctype), si->getRowLower(), si->getRowUpper(), probname);

	/** free memory */
	FREE_PTR(mat);
	FREE_ARRAY_PTR(ctype);
}

/** load problem */
void SolverInterfaceScip::loadProblem(
		CoinPackedMatrix * mat,
		const double * collb,
		const double * colub,
		const double * obj,
		const char * ctype,
		const double * rowlb,
		const double * rowub,
		const char * probname)
{
	assert(mat->getNumCols() > 0);
	char varname[128];
	char conname[128];
	/** create problem */
	SCIP_CALL_ABORT(SCIPcreateProbBasic(scip_, probname));

	/** allocate memory */
	vars_ = new SCIP_Var * [mat->getNumCols()];
	clbd_ = new double [mat->getNumCols()];
	cubd_ = new double [mat->getNumCols()];
	conss_ = new SCIP_CONS * [mat->getNumRows()];
	rlbd_ = new double [mat->getNumRows()];
	rubd_ = new double [mat->getNumRows()];

	/** add variables */
	for (int j = 0; j < mat->getNumCols(); ++j)
	{
		sprintf(varname, "var%d", j);
		SCIP_VAR * var = NULL;
		SCIP_VARTYPE vartype = ctype[j] == 'C' ? SCIP_VARTYPE_CONTINUOUS : SCIP_VARTYPE_INTEGER;
		SCIP_CALL_ABORT(SCIPcreateVar(scip_, &var, varname,
				collb[j],
				colub[j],
				obj[j],
				vartype,
				true,
				false,
				NULL, NULL, NULL, NULL, NULL));
		SCIP_CALL_ABORT(SCIPaddVar(scip_, var));
		vars_[j] = var;
		nvars_++;
		clbd_[j] = collb[j];
		cubd_[j] = colub[j];
	}

	/** add constraints */
	if (mat->isColOrdered())
		mat->reverseOrdering();
	for (int i = 0; i < mat->getNumRows(); ++i)
	{
		sprintf(conname, "con%d", i);
		int nvars = mat->getVectorSize(i);
		const int * indices = mat->getIndices() + mat->getVectorStarts()[i];
		const double * elements = mat->getElements() + mat->getVectorStarts()[i];

		SCIP_CONS * cons = NULL;
		SCIP_CALL_ABORT(SCIPcreateConsBasicLinear(scip_, &cons, conname, 0, NULL, NULL, rowlb[i], rowub[i]));

		for (int j = 0; j < nvars; ++j)
		{
			SCIP_CALL_ABORT(SCIPaddCoefLinear(scip_, cons, vars_[indices[j]], elements[j]));
		}

		/** add constraint */
		SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
		conss_[i] = cons;
		nconss_++;
		rlbd_[i] = rowlb[i];
		rubd_[i] = rowub[i];
	}
}

/** add row */
void SolverInterfaceScip::addRow(int size, const int * indices, const double * vals, double lb, double ub)
{
	assert(scip_);

	int nconss = SCIPgetNOrigConss(scip_);

	/** backup row info */
	SCIP_CONS ** conss = new SCIP_CONS * [nconss];
	double * rlbd = new double [nconss];
	double * rubd = new double [nconss];
	for (int i = 0; i < nconss; ++i)
	{
		conss[i] = conss_[i];
		conss_[i] = NULL;
		rlbd[i] = rlbd_[i];
		rubd[i] = rubd_[i];
	}
	FREE_ARRAY_PTR(conss_);
	FREE_ARRAY_PTR(rlbd_);
	FREE_ARRAY_PTR(rubd_);

	/** allocate with larger size */
	conss_ = new SCIP_CONS * [nconss+1];
	rlbd_ = new double [nconss+1];
	rubd_ = new double [nconss+1];

	/** restore row info */
	for (int i = 0; i < nconss; ++i)
	{
		conss_[i] = conss[i];
		conss[i] = NULL;
		rlbd_[i] = rlbd[i];
		rubd_[i] = rubd[i];
	}
	FREE_ARRAY_PTR(conss);
	FREE_ARRAY_PTR(rlbd);
	FREE_ARRAY_PTR(rubd);

	SCIP_CONS * cons = NULL;
	SCIP_CALL_ABORT(SCIPcreateConsBasicLinear(scip_, &cons, "cons", 0, NULL, NULL, lb, ub));

	/** add coefficients */
	for (int j = 0; j < size; ++j)
		SCIP_CALL_ABORT(SCIPaddCoefLinear(scip_, cons, vars_[indices[j]], vals[j]));

	/** add constraint */
	SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
	conss_[nconss] = cons;
	nconss_++;
	rlbd_[nconss] = lb;
	rubd_[nconss] = ub;
}

/** delete row */
void SolverInterfaceScip::delRow(int index)
{
	/** delete row */
	SCIP_CONS * cons = conss_[index];
	SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
	SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &cons));

	/** backup row info */
	SCIP_CONS ** conss = new SCIP_CONS * [nconss_];
	double * rlbd = new double [nconss_];
	double * rubd = new double [nconss_];
	CoinCopyN(conss_, nconss_, conss);
	CoinCopyN(rlbd_, nconss_, rlbd);
	CoinCopyN(rubd_, nconss_, rubd);
	FREE_ARRAY_PTR(conss_);
	FREE_ARRAY_PTR(rlbd_);
	FREE_ARRAY_PTR(rubd_);

	/** decrease number of rows */
	nconss_--;

	/** allocate with smaller size */
	conss_ = new SCIP_CONS * [nconss_];
	rlbd_ = new double [nconss_];
	rubd_ = new double [nconss_];
	if (index > 0)
	{
		CoinCopyN(conss, index, conss_);
		CoinCopyN(rlbd, index, rlbd_);
		CoinCopyN(rubd, index, rubd_);
	}
	if (index < nconss_)
	{
		CoinCopyN(conss + index + 1, nconss_ - index, conss_ + index);
		CoinCopyN(rlbd + index + 1, nconss_ - index, rlbd_ + index);
		CoinCopyN(rubd + index + 1, nconss_ - index, rubd_ + index);
	}
	FREE_ARRAY_PTR(conss);
	FREE_ARRAY_PTR(rlbd);
	FREE_ARRAY_PTR(rubd);
}

/** add constraint handler */
void SolverInterfaceScip::addConstraintHandler(
		scip::ObjConshdlr * objconshdlr,
		bool deleteobject,
		bool isDual)
{
	assert(scip_);

	/** include constraint handler */
	SCIP_CALL_ABORT(SCIPincludeObjConshdlr(scip_, objconshdlr, deleteobject));

	/* create constraint */
	SCIP_CONS * cons = NULL;
	if (!isDual)
	{
		SCIP_CALL_ABORT(SCIPcreateConsBenders(scip_, &cons, "Benders"));
		SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
		SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &cons));
	}
	else
	{
		SCIP_CALL_ABORT(SCIPcreateConsBendersDd(scip_, &cons, "BendersDd"));
		SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
		SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &cons));
	}
}

/** find constriant handler */
scip::ObjConshdlr * SolverInterfaceScip::findObjConshdlr(const char * name)
{
	return SCIPfindObjConshdlr(scip_, name);
}

/** add branch rule */
void SolverInterfaceScip::addBranchrule(
		scip::ObjBranchrule * objbranchrule,
		bool deleteobject)
{
	assert(scip_);
	/** include constraint handler */
	SCIP_CALL_ABORT(SCIPincludeObjBranchrule(scip_, objbranchrule, deleteobject));
}

/** find branch rule */
scip::ObjBranchrule * SolverInterfaceScip::findObjBranchrule(const char * name)
{
	return SCIPfindObjBranchrule(scip_, name);
}

#if 0
/** enable upper bounding cuts */
void SolverInterfaceScip::enableUpperBoundingCuts(bool yesNo, , bool hasUb)
{
	SCIP_CONS * cons = SCIPfindCons(scip_, "BendersDd");
	if (yesNo)
	{
		if (cons == NULL)
		{
			SCIP_CALL_ABORT(SCIPcreateConsBendersDd(scip_, &cons, "BendersDd", hasUb));
			SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
			SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &cons));
		}
	}
	else if (cons)
	{
//		printf("Disable upper bounding cuts %p", cons);
		SCIP_CALL_ABORT(SCIPdelCons(scip_, cons));
//		printf(" %p\n", cons);
	}
}
#endif

/** solve */
void SolverInterfaceScip::solve()
{
	/** solve */
	SCIP_CALL_ABORT(SCIPsolve(scip_));
#if DSP_DEBUG
	int numDdCuts = 0;
	int numPoolCuts = SCIPgetNPoolCuts(scip_);
	SCIP_CUT ** poolcuts = SCIPgetPoolCuts(scip_);
	for (int i = numPoolCuts - 1; i >= 0; --i)
	{
		/** retrieve row */
		SCIP_ROW * poolcutrow = SCIPcutGetRow(poolcuts[i]);

		/** benders? */
		if (strcmp(SCIProwGetName(poolcutrow), "bendersDd") == 0)
			numDdCuts++;
	}
	DSPdebugMessage("Number of pool cuts %d (DdCuts %d)\n", numPoolCuts, numDdCuts);
#endif

	/** get best solution */
	SCIP_SOL * sol = SCIPgetBestSol(scip_);
	if (sol)
	{
		/** allocate memory */
		FREE_ARRAY_PTR(solution_);
		solution_ = new double [SCIPgetNOrigVars(scip_)];

		/** get solution values */
		SCIP_CALL_ABORT(SCIPgetSolVals(scip_, sol, SCIPgetNOrigVars(scip_), vars_, solution_));
	}
}

/** get solution values */
const double * SolverInterfaceScip::getSolution()
{
	return solution_;
}

/** get array of solutions; returns the number of solutions */
int SolverInterfaceScip::getSolutions(double *** solutions)
{
	assert(*solutions == NULL);

	/** get number of solutions */
	int nsols = SCIPgetNSols(scip_);

	if (nsols == 0) return 0;

	/** get solutions */
	SCIP_SOL ** sols = SCIPgetSols(scip_);

	/** allocate memory */
	*solutions = new double * [nsols];
	for (int i = 0; i < nsols; ++i)
	{
		/** allocate memory */
		(*solutions)[i] = new double [SCIPgetNOrigVars(scip_)];

		/** get solution values */
		SCIP_CALL_ABORT(SCIPgetSolVals(scip_, sols[i], SCIPgetNOrigVars(scip_), vars_, (*solutions)[i]));
	}

	return nsols;
}

/** get cuts */
const OsiCuts * SolverInterfaceScip::getCuts()
{
	scip::ObjConshdlr * conshdlr = SCIPfindObjConshdlr(scip_, "BendersDd");
	if (!conshdlr) return NULL;

	SCIPconshdlrBendersDd * conshdlrDd = dynamic_cast<SCIPconshdlrBendersDd*>(conshdlr);
	if (!conshdlrDd) return NULL;

	return conshdlrDd->getCutsAdded();
}

/** set cuts */
void SolverInterfaceScip::setCuts(OsiCuts * cuts)
{
	scip::ObjConshdlr * conshdlr = SCIPfindObjConshdlr(scip_, "BendersDd");
	if (!conshdlr) return;

	SCIPconshdlrBendersDd * conshdlrDd = dynamic_cast<SCIPconshdlrBendersDd*>(conshdlr);
	if (!conshdlrDd) return;

	conshdlrDd->setCutsToAdd(cuts);
}

/** clear cuts */
void SolverInterfaceScip::clearCuts()
{
	scip::ObjConshdlr * conshdlr = SCIPfindObjConshdlr(scip_, "BendersDd");
	if (!conshdlr) return;

	SCIPconshdlrBendersDd * conshdlrDd = dynamic_cast<SCIPconshdlrBendersDd*>(conshdlr);
	if (!conshdlrDd) return;

	conshdlrDd->clearCutsAdded();
}

/** solution status */
STO_RTN_CODE SolverInterfaceScip::getStatus()
{
	/** solution status */
	switch (SCIPgetStatus(scip_))
	{
	case SCIP_STATUS_OPTIMAL:
		return STO_STAT_OPTIMAL;
	case SCIP_STATUS_INFEASIBLE:
		return STO_STAT_PRIM_INFEASIBLE;
	case SCIP_STATUS_UNBOUNDED:
		return STO_STAT_DUAL_INFEASIBLE;
	case SCIP_STATUS_NODELIMIT:
		return STO_STAT_STOPPED_NODE;
	case SCIP_STATUS_TIMELIMIT:
		return STO_STAT_STOPPED_TIME;
	case SCIP_STATUS_GAPLIMIT:
		return STO_STAT_STOPPED_GAP;
	default:
		printf("Unexpected solution status %d\n", SCIPgetStatus(scip_));
		writeMps("unexpectedStatus.mps");	
		return STO_STAT_UNKNOWN;
	}
}

/** transform problem */
void SolverInterfaceScip::transformProb()
{
	SCIP_CALL_ABORT(SCIPtransformProb(scip_));
}

/** transform and presolve */
void SolverInterfaceScip::presolve()
{
	SCIP_CALL_ABORT(SCIPpresolve(scip_));
}

/** free solve */
void SolverInterfaceScip::freeSolve(bool restart)
{
	if (SCIPgetStage(scip_) == SCIP_STAGE_SOLVED ||
		SCIPgetStage(scip_) == SCIP_STAGE_SOLVING ||
		SCIPgetStage(scip_) == SCIP_STAGE_TRANSFORMED)
		SCIP_CALL_ABORT(SCIPfreeSolve(scip_, restart));
}

/** free trensformed problem */
void SolverInterfaceScip::freeTransform()
{
	if (SCIPgetStage(scip_) == SCIP_STAGE_INIT ||
			SCIPgetStage(scip_) == SCIP_STAGE_PROBLEM ||
			SCIPgetStage(scip_) == SCIP_STAGE_TRANSFORMED ||
			SCIPgetStage(scip_) == SCIP_STAGE_PRESOLVING ||
			SCIPgetStage(scip_) == SCIP_STAGE_PRESOLVED ||
			SCIPgetStage(scip_) == SCIP_STAGE_SOLVING ||
			SCIPgetStage(scip_) == SCIP_STAGE_SOLVED)
		SCIP_CALL_ABORT(SCIPfreeTransform(scip_));
}

/** get print out level */
int SolverInterfaceScip::getPrintLevel()
{
	assert(scip_);
	int level;
	SCIP_CALL_ABORT(SCIPgetIntParam(scip_, "display/verblevel", &level));
	return level;
}

/** set solution */
void SolverInterfaceScip::setSolution(double * solution)
{
	if (solution == NULL) return;

	if (SCIPgetStage(scip_) == SCIP_STAGE_PROBLEM)
		SCIP_CALL_ABORT(SCIPtransformProb(scip_));

	SCIP_Sol * sol = NULL;
	SCIP_Bool stored;

	/** create solution */
	SCIP_CALL_ABORT(SCIPcreateSol(scip_, &sol, NULL));

	/** set solution values */
	SCIP_CALL_ABORT(SCIPsetSolVals(scip_, sol, SCIPgetNOrigVars(scip_), vars_, solution));

	/** check solution and free if infeasible */
	SCIP_CALL_ABORT(SCIPtrySolFree(scip_, &sol, false, true, true, true, &stored));
	//printf("DEBUG: Is solution added? %s\n", stored ? "yes" : "no");
}

/** set objective coefficients */
void SolverInterfaceScip::setObjCoef(double * obj)
{
	/** free transform */
	freeTransform();
	int nvars = SCIPgetNOrigVars(scip_);
	for (int j = 0; j < nvars; ++j)
		SCIP_CALL_ABORT(SCIPchgVarObj(scip_, vars_[j], obj[j]));
}

/** set objective coefficients */
void SolverInterfaceScip::setObjSense(int sense)
{
	assert(scip_);
	SCIP_CALL_ABORT(SCIPsetObjsense(scip_, sense > 0 ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE));
}

/** set column bounds */
void SolverInterfaceScip::setColLower(int index, double lb)
{
	/** free transform */
	freeTransform();
	SCIP_CALL_ABORT(SCIPchgVarLb(scip_, vars_[index], lb));
	clbd_[index] = lb;
}

/** set column bounds */
void SolverInterfaceScip::setColUpper(int index, double ub)
{
	/** free transform */
	freeTransform();
	SCIP_CALL_ABORT(SCIPchgVarUb(scip_, vars_[index], ub));
	cubd_[index] = ub;
}

/** set column bounds */
void SolverInterfaceScip::setColBounds(int index, double lb, double ub)
{
	setColLower(index, lb);
	setColUpper(index, ub);
}

/** set row bounds */
void SolverInterfaceScip::setRowLower(int index, double lb)
{
	/** free transform */
	freeTransform();
	assert(conss_);
	assert(rlbd_);
	SCIP_CALL_ABORT(SCIPchgLhsLinear(scip_, conss_[index], lb));
	rlbd_[index] = lb;
}

/** set row bounds */
void SolverInterfaceScip::setRowUpper(int index, double ub)
{
	/** free transform */
	freeTransform();
	assert(conss_);
	assert(rubd_);
	SCIP_CALL_ABORT(SCIPchgRhsLinear(scip_, conss_[index], ub));
	rubd_[index] = ub;
}

/** set row bounds */
void SolverInterfaceScip::setRowBounds(int index, double lb, double ub)
{
	setRowLower(index, lb);
	setRowUpper(index, ub);
}

/** set print out level */
void SolverInterfaceScip::setPrintLevel(int level)
{
	assert(scip_);
	SCIP_CALL_ABORT(SCIPsetIntParam(scip_, "display/verblevel", level));
}

/** set node limit */
void SolverInterfaceScip::setNodeLimit(int limit)
{
	SCIP_CALL_ABORT(SCIPsetLongintParam(scip_, "limits/nodes", limit));
}

/** set clock type */
void SolverInterfaceScip::setClockType(int type)
{
	SCIP_CALL_ABORT(SCIPsetIntParam(scip_, "timing/clocktype", type));
}

void SolverInterfaceScip::setTimeLimit(double sec)
{
	SCIP_CALL_ABORT(SCIPsetRealParam(scip_, "limits/time", sec));
}

/** dual reduction */
SCIP_RETCODE SolverInterfaceScip::setDualReduction(bool yesNo)
{
	SCIP_CALL(SCIPsetBoolParam(scip_, "constraints/indicator/dualreductions", yesNo));
	SCIP_CALL(SCIPsetBoolParam(scip_, "constraints/setppc/dualpresolving", yesNo));
	SCIP_CALL(SCIPsetBoolParam(scip_, "constraints/linear/dualpresolving", yesNo));
	SCIP_CALL(SCIPsetBoolParam(scip_, "constraints/abspower/dualpresolve", yesNo));
	SCIP_CALL(SCIPsetBoolParam(scip_, "constraints/and/dualpresolving", yesNo));
	SCIP_CALL(SCIPsetBoolParam(scip_, "constraints/cumulative/dualpresolve", yesNo));
	return SCIP_OKAY;
}

/** write MPS file */
void SolverInterfaceScip::writeMps(const char * filename)
{
	assert(scip_);
	SCIP_CALL_ABORT(SCIPwriteOrigProblem(scip_, filename, "mps", false));

}

