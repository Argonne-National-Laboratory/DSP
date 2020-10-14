/*
 * OsiScipSolverInterface.cpp
 *
 *  Created on: Sep 9, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include <memory>

#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"

#include "scip/scipdefplugins.h"

#include "Utility/DspMessage.h"
#include "SolverInterface/OsiScipSolverInterface.hpp"

void OsiScipSolverInterface::writeMps(const char* filename,
		const char* extension, double objSense) const {
	assert(scip_);
	SCIP_CALL_ABORT(SCIPwriteOrigProblem(scip_, filename, extension, false));
}

void OsiScipSolverInterface::initialSolve() {
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
		/** get solution values */
		solution_.resize(SCIPgetNOrigVars(scip_));
		SCIP_CALL_ABORT(SCIPgetSolVals(scip_, sol, SCIPgetNOrigVars(scip_), &vars_[0], &solution_[0]));
	}
}

void OsiScipSolverInterface::resolve() {
	initialSolve();
}

void OsiScipSolverInterface::branchAndBound() {
	initialSolve();
}

bool OsiScipSolverInterface::isAbandoned() const {
	return (SCIPgetStatus(scip_) == SCIP_STATUS_UNKNOWN);
}

bool OsiScipSolverInterface::isProvenOptimal() const {
	return (SCIPgetStatus(scip_) == SCIP_STATUS_OPTIMAL);
}

bool OsiScipSolverInterface::isProvenPrimalInfeasible() const {
	return (SCIPgetStatus(scip_) == SCIP_STATUS_INFEASIBLE);
}

bool OsiScipSolverInterface::isProvenDualInfeasible() const {
	return (SCIPgetStatus(scip_) == SCIP_STATUS_UNBOUNDED);
}

bool OsiScipSolverInterface::isIterationLimitReached() const {
	return (SCIPgetStatus(scip_) == SCIP_STATUS_NODELIMIT);
}

CoinWarmStart* OsiScipSolverInterface::getEmptyWarmStart() const {
	CoinError("Not supported.", "getEmptyWarmStart", "OsiScipSolverInterface");
	return NULL;
}

CoinWarmStart* OsiScipSolverInterface::getWarmStart() const {
	CoinError("Not supported.", "getWarmStart", "OsiScipSolverInterface");
	return NULL;
}

bool OsiScipSolverInterface::setWarmStart(const CoinWarmStart* warmstart) {
	CoinError("Not supported.", "setWarmStart", "OsiScipSolverInterface");
	return false;
}

int OsiScipSolverInterface::getNumCols() const {
	return SCIPgetNOrigVars(scip_);
}

int OsiScipSolverInterface::getNumRows() const {
	return SCIPgetNOrigConss(scip_);
}

int OsiScipSolverInterface::getNumElements() const {
	return mat_->getNumElements();
}

const double* OsiScipSolverInterface::getColLower() const {
	return clbd_.data();
}

const double* OsiScipSolverInterface::getColUpper() const {
	return cubd_.data();
}

const char* OsiScipSolverInterface::getRowSense() const {
	return sense_.data();
}

const double* OsiScipSolverInterface::getRightHandSide() const {
	return rhs_.data();
}

const double* OsiScipSolverInterface::getRowRange() const {
	return range_.data();
}

const double* OsiScipSolverInterface::getRowLower() const {
	return rlbd_.data();
}

const double* OsiScipSolverInterface::getRowUpper() const {
	return rubd_.data();
}

const double* OsiScipSolverInterface::getObjCoefficients() const {
	return obj_.data();
}

double OsiScipSolverInterface::getObjSense() const {
	return SCIPgetObjsense(scip_) == SCIP_OBJSENSE_MINIMIZE ? 1 : -11;
}

bool OsiScipSolverInterface::isContinuous(int colIndex) const {
	return (ctype_[colIndex] == 'C');
}

const CoinPackedMatrix* OsiScipSolverInterface::getMatrixByRow() const {
	if (mat_->isColOrdered())
		mat_->reverseOrdering();
	return mat_;
}

const CoinPackedMatrix* OsiScipSolverInterface::getMatrixByCol() const {
	if (mat_->isColOrdered() == false)
		mat_->reverseOrdering();
	return mat_;
}

double OsiScipSolverInterface::getInfinity() const {
	return SCIPinfinity(scip_);
}

const double* OsiScipSolverInterface::getColSolution() const {
	return solution_.data();
}

const double* OsiScipSolverInterface::getRowPrice() const {
	return price_.data();
}

const double* OsiScipSolverInterface::getReducedCost() const {
	return reduced_.data();
}

const double* OsiScipSolverInterface::getRowActivity() const {
	return activity_.data();
}

double OsiScipSolverInterface::getObjValue() const {
	return SCIPgetPrimalbound(scip_);
}

int OsiScipSolverInterface::getIterationCount() const {
	return SCIPgetNLPIterations(scip_);
}
//------------------------------------------------------------------
int OsiScipSolverInterface::getNumNodes() const
{
	return SCIPgetNNodes(scip_);
}
//------------------------------------------------------------------
double OsiScipSolverInterface::getBestDualBound() const
{
  return SCIPgetDualbound(scip_);
}

std::vector<double*> OsiScipSolverInterface::getDualRays(
		int maxNumRays,
		bool fullRay) const {
	/** TODO: Implement */
	CoinError("Not implemented.", "getDualRays", "OsiScipSolverInterface");
	return std::vector<double*>();
}

std::vector<double*> OsiScipSolverInterface::getPrimalRays(int maxNumRays) const {
	/** TODO: Implement */
	CoinError("Not implemented.", "getPrimalRays", "OsiScipSolverInterface");
	return std::vector<double*>();
}

void OsiScipSolverInterface::setObjCoeff(
		int elementIndex,
		double elementValue) {
	freeTransform();
	SCIP_CALL_ABORT(SCIPchgVarObj(scip_, vars_[elementIndex], elementValue));
	obj_[elementIndex] = elementValue;
}

void OsiScipSolverInterface::setObjSense(double s) {
	assert(scip_);
	SCIP_CALL_ABORT(SCIPsetObjsense(scip_, s > 0 ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE));
}

void OsiScipSolverInterface::setColLower(int elementIndex,
		double elementValue) {
	freeTransform();
	SCIP_CALL_ABORT(SCIPchgVarLb(scip_, vars_[elementIndex], elementValue));
	clbd_[elementIndex] = elementValue;
}

void OsiScipSolverInterface::setColUpper(int elementIndex,
		double elementValue) {
	freeTransform();
	SCIP_CALL_ABORT(SCIPchgVarUb(scip_, vars_[elementIndex], elementValue));
	clbd_[elementIndex] = elementValue;
}

void OsiScipSolverInterface::setRowLower(int elementIndex,
		double elementValue) {
	freeTransform();
	assert(conss_[elementIndex]);
	SCIP_CALL_ABORT(SCIPchgLhsLinear(scip_, conss_[elementIndex], elementValue));

	rlbd_[elementIndex] = elementValue;
	convertBoundToSense(rlbd_[elementIndex], rubd_[elementIndex], 
		sense_[elementIndex], rhs_[elementIndex], range_[elementIndex]);
}

void OsiScipSolverInterface::setRowUpper(int elementIndex,
		double elementValue) {
	freeTransform();
	assert(conss_[elementIndex]);
	SCIP_CALL_ABORT(SCIPchgRhsLinear(scip_, conss_[elementIndex], elementValue));

	rubd_[elementIndex] = elementValue;
	convertBoundToSense(rlbd_[elementIndex], rubd_[elementIndex], 
		sense_[elementIndex], rhs_[elementIndex], range_[elementIndex]);
}

void OsiScipSolverInterface::setRowType(int index, char sense,
		double rightHandSide, double range) {
	sense_[index] = sense;
	range_[index] = range;
	rhs_[index] = rightHandSide;
	convertSenseToBound(sense_[index], rhs_[index], range_[index], 
		rlbd_[index], rubd_[index]);
}

void OsiScipSolverInterface::setColSolution(const double* colsol) {
	if (solution_.size() != vars_.size()) 
		return;

	if (SCIPgetStage(scip_) == SCIP_STAGE_PROBLEM)
		SCIP_CALL_ABORT(SCIPtransformProb(scip_));

	SCIP_Sol * sol = NULL;
	SCIP_Bool stored;

	/** create solution */
	SCIP_CALL_ABORT(SCIPcreateSol(scip_, &sol, NULL));

	/** set solution values */
	SCIP_CALL_ABORT(SCIPsetSolVals(scip_, sol, SCIPgetNOrigVars(scip_), &vars_[0], &solution_[0]));

	/** check solution and free if infeasible */
#if SCIP_VERSION_MAJOR < 6
	SCIP_CALL_ABORT(SCIPtrySolFree(scip_, &sol, true, true, true, true, &stored));
#else
	SCIP_CALL_ABORT(SCIPtrySolFree(scip_, &sol, true, true, true, true, true, &stored));
#endif
}

void OsiScipSolverInterface::setRowPrice(const double* rowprice) {
	/** TODO: Not supported */
	CoinError("Not supported.", "setRowPrice", "OsiScipSolverInterface");
}

//-----------------------------------------------------------------------------
void OsiScipSolverInterface::setTimeLimit(double t)
{
  SCIPsetRealParam(scip_, "limits/time", t);
}
//-----------------------------------------------------------------------------
void OsiScipSolverInterface::setNodeLimit(int n)
{
  SCIPsetLongintParam(scip_, "limits/nodes", n);
}
//-----------------------------------------------------------------------------
void OsiScipSolverInterface::setMipRelGap(double gap)
{
  SCIPsetRealParam(scip_, "limits/gap", gap);
}

void OsiScipSolverInterface::setContinuous(int index) {
	if (ctype_[index] != 'C')
	{
		ctype_[index] = 'C';
		SCIP_Bool infeasible;
		SCIP_CALL_ABORT(SCIPchgVarType(scip_, vars_[index], SCIP_VARTYPE_CONTINUOUS, &infeasible));
		numberIntegers_--;
	}
}

void OsiScipSolverInterface::setInteger(int index) {                    
	if (ctype_[index] == 'C')
	{
		SCIP_VARTYPE scip_type;
		SCIP_Bool infeasible;
		if (clbd_[index] == 0.0 && cubd_[index] == 1.0)
		{
			ctype_[index] = 'B';
			scip_type = SCIP_VARTYPE_BINARY;
		}
		else
		{
			ctype_[index] = 'I';
			scip_type = SCIP_VARTYPE_INTEGER;
		}
		SCIP_CALL_ABORT(SCIPchgVarType(scip_, vars_[index], scip_type, &infeasible));
		numberIntegers_++;
	}
}

void OsiScipSolverInterface::addCol(
		const CoinPackedVectorBase& vec,
		const double collb,
		const double colub,
		const double obj) {
	assert(scip_);

	/** create new variable */
	SCIP_VAR * var;
	SCIP_CALL_ABORT(SCIPcreateVar(scip_, &var, "var", collb, colub, obj,
								  SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL));

	/** add the variable */
	SCIP_CALL_ABORT(SCIPaddVar(scip_, var));

	/** add coefficients */
	for (int i = 0; i < vec.getNumElements(); ++i) {
		SCIP_CALL_ABORT(SCIPaddCoefLinear(scip_, conss_[vec.getIndices()[i]], var, vec.getElements()[i]));
	}

	/** add column */
	vars_.push_back(var);
	clbd_.push_back(collb);
	cubd_.push_back(colub);
	ctype_.push_back('C');
	obj_.push_back(obj);
	nvars_++;
}

void OsiScipSolverInterface::deleteCols(const int num, const int* colIndices) {

	SCIP_Bool deleted;
	std::vector<int> inds;
	inds.resize(num);

	freeTransform();

	for (int j = 0; j < num; ++j) 
	{
		SCIP_CALL_ABORT(SCIPdelVar(scip_, vars_[colIndices[j]], &deleted));
		//SCIP_CALL_ABORT(SCIPreleaseVar(scip_, &vars_[colIndices[j]]));
		inds[j] = colIndices[j];
	}

	std::sort(inds.begin(), inds.end());
	std::reverse(inds.begin(), inds.end());
	for (int j = 0; j < num; ++j)
	{
		vars_.erase(vars_.begin() + inds[j]);
		clbd_.erase(clbd_.begin() + inds[j]);
		cubd_.erase(cubd_.begin() + inds[j]);
		ctype_.erase(ctype_.begin() + inds[j]);
	}
}

void OsiScipSolverInterface::addRow(
		const CoinPackedVectorBase& vec,
		const double rowlb,
		const double rowub) {
	assert(scip_);

	SCIP_CONS * cons = NULL;
	SCIP_CALL_ABORT(SCIPcreateConsBasicLinear(scip_, &cons, "cons", 0, NULL, NULL, rowlb, rowub));

	/** add coefficients */
	for (int j = 0; j < vec.getNumElements(); ++j)
		SCIP_CALL_ABORT(SCIPaddCoefLinear(scip_, cons, vars_[vec.getIndices()[j]], vec.getElements()[j]));

	/** add constraint */
	SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
	conss_.push_back(cons);
	rlbd_.push_back(rowlb);
	rubd_.push_back(rowub);

	char sense; double right, range;
	convertBoundToSense(rowlb, rowub, sense, right, range);
	sense_.push_back(sense);
	rhs_.push_back(right);
	range_.push_back(range);
	nconss_++;
}

void OsiScipSolverInterface::addRow(
		const CoinPackedVectorBase& vec,
		const char rowsen,
		const double rowrhs,
		const double rowrng) {
	double rowlb, rowub;
	convertSenseToBound(rowsen, rowrhs, rowrng, rowlb, rowub);
	addRow(vec, rowlb, rowub);
}

void OsiScipSolverInterface::deleteRows(const int num, const int* rowIndices) {
	std::vector<int> inds;
	inds.resize(num);

	freeTransform();

	for (int j = 0; j < num; ++j) 
	{
		SCIP_CALL_ABORT(SCIPdelCons(scip_, conss_[rowIndices[j]]));
		SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &conss_[rowIndices[j]]));
		inds[j] = rowIndices[j];
	}

	std::sort(inds.begin(), inds.end());
	std::reverse(inds.begin(), inds.end());
	for (int j = 0; j < num; ++j)
	{
		conss_.erase(conss_.begin() + inds[j]);
		rlbd_.erase(rlbd_.begin() + inds[j]);
		rubd_.erase(rubd_.begin() + inds[j]);
		sense_.erase(sense_.begin() + inds[j]);
		rhs_.erase(rhs_.begin() + inds[j]);
		range_.erase(range_.begin() + inds[j]);
	}
}

void OsiScipSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
		const double* collb, const double* colub, const double* obj,
		const double* rowlb, const double* rowub) {
	assert(matrix.getNumCols() > 0);
	char varname[128];
	char conname[128];

	/** create problem */
	SCIP_CALL_ABORT(SCIPcreateProbBasic(scip_, "Dsp"));

	/** allocate memory */
	if (mat_) delete mat_;
	mat_ = new CoinPackedMatrix(matrix);
	vars_.resize(mat_->getNumCols());
	obj_.resize(mat_->getNumCols());
	clbd_.resize(mat_->getNumCols());
	cubd_.resize(mat_->getNumCols());
	conss_.resize(mat_->getNumRows());
	rlbd_.resize(mat_->getNumRows());
	rubd_.resize(mat_->getNumRows());
	sense_.resize(mat_->getNumRows());
	rhs_.resize(mat_->getNumRows());
	range_.resize(mat_->getNumRows());
	ctype_.resize(mat_->getNumCols(), 'C');
	nvars_ = mat_->getNumCols();

	/** add variables */
	for (int j = 0; j < mat_->getNumCols(); ++j)
	{
		sprintf(varname, "var%d", j);
		SCIP_VAR * var = NULL;
		SCIP_VARTYPE vartype = SCIP_VARTYPE_CONTINUOUS;
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
		obj_[j] = obj[j];
		clbd_[j] = collb[j];
		cubd_[j] = colub[j];
	}

	/** add constraints */
	if (mat_->isColOrdered())
		mat_->reverseOrdering();
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		sprintf(conname, "con%d", i);
		int nvars = mat_->getVectorSize(i);
		const int * indices = mat_->getIndices() + mat_->getVectorStarts()[i];
		const double * elements = mat_->getElements() + mat_->getVectorStarts()[i];

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
		convertBoundToSense(rlbd_[i], rubd_[i], sense_[i], rhs_[i], range_[i]);
	}
}

void OsiScipSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
		double*& collb, double*& colub, double*& obj, double*& rowlb,
		double*& rowub) {
	loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
}

void OsiScipSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
		const double* collb, const double* colub, const double* obj,
		const char* rowsen, const double* rowrhs, const double* rowrng) {
	assert(matrix.getNumCols() > 0);
	char varname[128];
	char conname[128];

	/** create problem */
	SCIP_CALL_ABORT(SCIPcreateProbBasic(scip_, "Dsp"));

	/** allocate memory */
	if (mat_) delete mat_;
	mat_ = new CoinPackedMatrix(matrix);
	vars_.resize(mat_->getNumCols());
	obj_.resize(mat_->getNumCols());
	clbd_.resize(mat_->getNumCols());
	cubd_.resize(mat_->getNumCols());
	conss_.resize(mat_->getNumRows());
	rlbd_.resize(mat_->getNumRows());
	rubd_.resize(mat_->getNumRows());
	sense_.resize(mat_->getNumRows());
	rhs_.resize(mat_->getNumRows());
	range_.resize(mat_->getNumRows());
	ctype_.resize(mat_->getNumCols(), 'C');
	nvars_ = mat_->getNumCols();

	/** add variables */
	for (int j = 0; j < mat_->getNumCols(); ++j)
	{
		sprintf(varname, "var%d", j);
		SCIP_VAR * var = NULL;
		SCIP_VARTYPE vartype = SCIP_VARTYPE_CONTINUOUS;
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
		obj_[j] = obj[j];
		clbd_[j] = collb[j];
		cubd_[j] = colub[j];
	}

	/** add constraints */
	if (mat_->isColOrdered())
		mat_->reverseOrdering();
	for (int i = 0; i < mat_->getNumRows(); ++i)
	{
		sprintf(conname, "con%d", i);
		int nvars = mat_->getVectorSize(i);
		const int * indices = mat_->getIndices() + mat_->getVectorStarts()[i];
		const double * elements = mat_->getElements() + mat_->getVectorStarts()[i];

		sense_[i] = rowsen[i];
		rhs_[i] = rowrhs[i];
		range_[i] = rowrng[i];
		convertSenseToBound(rlbd_[i], rubd_[i], sense_[i], rhs_[i], range_[i]);

		SCIP_CONS * cons = NULL;
		SCIP_CALL_ABORT(SCIPcreateConsBasicLinear(scip_, &cons, conname, 0, NULL, NULL, rlbd_[i], rubd_[i]));

		for (int j = 0; j < nvars; ++j)
		{
			SCIP_CALL_ABORT(SCIPaddCoefLinear(scip_, cons, vars_[indices[j]], elements[j]));
		}

		/** add constraint */
		SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
		conss_[i] = cons;
		nconss_++;
	}
}

void OsiScipSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
		double*& collb, double*& colub, double*& obj, char*& rowsen,
		double*& rowrhs, double*& rowrng) {
	loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
}

void OsiScipSolverInterface::loadProblem(const int numcols, const int numrows,
		const CoinBigIndex* start, const int* index, const double* value,
		const double* collb, const double* colub, const double* obj,
		const double* rowlb, const double* rowub) {
	/** column length */
	std::vector<int> len(numcols);
	for (int j = 0; j < numcols; ++j)
		len[j] = start[j+1] - start[j];
	/** create matrix */
	CoinPackedMatrix* matrix = new CoinPackedMatrix(true, numrows, numcols, start[numcols], value, index, start, len.data());
	/** load problem */
	loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
	/** release memory */
	delete matrix;
	matrix = NULL;
}

void OsiScipSolverInterface::loadProblem(const int numcols, const int numrows,
		const CoinBigIndex* start, const int* index, const double* value,
		const double* collb, const double* colub, const double* obj,
		const char* rowsen, const double* rowrhs, const double* rowrng) {
	/** column length */
	std::vector<int> len(numcols);
	for (int j = 0; j < numcols; ++j)
		len[j] = start[j+1] - start[j];
	/** create matrix */
	CoinPackedMatrix* matrix = new CoinPackedMatrix(true, numrows, numcols, start[numcols], value, index, start, len.data());
	/** load problem */
	loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
	/** release memory */
	delete matrix;
	matrix = NULL;
}

void OsiScipSolverInterface::initialize()
{
	if (scip_ != NULL)
	{
		SCIP_CALL_ABORT(SCIPfree(&scip_));
		scip_ = NULL;
	}
	SCIP_CALL_ABORT(SCIPcreate(&scip_));
	SCIP_CALL_ABORT(SCIPincludeDefaultPlugins(scip_));
}

void OsiScipSolverInterface::finalize()
{
	if (scip_)
	{
		for (int j = 0; j < nvars_; ++j)
			SCIP_CALL_ABORT(SCIPreleaseVar(scip_, &vars_[j]));
		for (int i = 0; i < nconss_; ++i)
			SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &conss_[i]));
		SCIP_CALL_ABORT(SCIPfree(&scip_));
		scip_ = NULL;
	}

	solution_.clear();
	price_.clear();
	reduced_.clear();
	activity_.clear();

	if (mat_) delete mat_;

	nvars_ = 0;
	vars_.clear();
	obj_.clear();
	clbd_.clear();
	cubd_.clear();
	ctype_.clear();

	nconss_ = 0;
	conss_.clear();
	rlbd_.clear();
	rubd_.clear();
	sense_.clear();
	rhs_.clear();
	range_.clear();
}

/** free trensformed problem */
void OsiScipSolverInterface::freeTransform()
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

// void OsiScipSolverInterface::addConstraintHandler(
// 	scip::ObjConshdlr * objconshdlr,
// 	bool deleteobject,
// 	bool isDual = false)
// {
// 	assert(scip_);

// 	/** include constraint handler */
// 	SCIP_CALL_ABORT(SCIPincludeObjConshdlr(scip_, objconshdlr, deleteobject));

// 	/* create constraint */
// 	SCIP_CONS * cons = NULL;
// 	if (!isDual)
// 	{
// 		SCIP_CALL_ABORT(SCIPcreateConsBenders(scip_, &cons, "Benders"));
// 		SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
// 		SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &cons));
// 	}
// 	else
// 	{
// 		SCIP_CALL_ABORT(SCIPcreateConsBenders(scip_, &cons, "BendersDd"));
// 		SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
// 		SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &cons));
// 	}
// }

/** find constriant handler */
scip::ObjConshdlr * OsiScipSolverInterface::findObjConshdlr(const char * name)
{
	return SCIPfindObjConshdlr(scip_, name);
}

/** add branch rule */
void OsiScipSolverInterface::addBranchrule(
		scip::ObjBranchrule * objbranchrule,
		bool deleteobject)
{
	assert(scip_);
	/** include constraint handler */
	SCIP_CALL_ABORT(SCIPincludeObjBranchrule(scip_, objbranchrule, deleteobject));
}

/** find branch rule */
scip::ObjBranchrule * OsiScipSolverInterface::findObjBranchrule(const char * name)
{
	return SCIPfindObjBranchrule(scip_, name);
}

OsiScipSolverInterface::OsiScipSolverInterface() :
scip_(NULL),
mat_(NULL),
nvars_(0),
nconss_(0) {
	initialize();
}

OsiScipSolverInterface* OsiScipSolverInterface::clone(bool copyData) const {
	return (new OsiScipSolverInterface(*this));
}

OsiScipSolverInterface::OsiScipSolverInterface(const OsiScipSolverInterface& rhs) :
OsiSolverInterface(rhs) {
	scip_ = NULL;
	mat_ = NULL;
	// Initialize SCIP
	initialize();
	/** load problem */
	loadProblem(*(rhs.mat_), &rhs.clbd_[0], &rhs.cubd_[0], &rhs.obj_[0], &rhs.rlbd_[0], &rhs.rubd_[0]);
	ctype_ = rhs.ctype_;
	/** copy solution data */
	solution_ = rhs.solution_;
	reduced_ = rhs.reduced_;
	price_ = rhs.price_;
	activity_ = rhs.activity_;
}

OsiScipSolverInterface& OsiScipSolverInterface::operator =(
		const OsiScipSolverInterface& rhs) {
	if (this != &rhs) {
		// Initialize SCIP
		scip_ = NULL;
		initialize();
		/** load problem */
		loadProblem(*(rhs.mat_), &rhs.clbd_[0], &rhs.cubd_[0], &rhs.obj_[0], &rhs.rlbd_[0], &rhs.rubd_[0]);
		ctype_ = rhs.ctype_;
		/** copy solution data */
		solution_ = rhs.solution_;
		reduced_ = rhs.reduced_;
		price_ = rhs.price_;
		activity_ = rhs.activity_;
	}
	return *this;
}

OsiScipSolverInterface::~OsiScipSolverInterface() {
	finalize();
}

void OsiScipSolverInterface::applyRowCut(const OsiRowCut& rc) {
	addRow(rc.row(), rc.lb(), rc.ub());
}

void OsiScipSolverInterface::applyColCut(const OsiColCut& cc) {
	CoinPackedVector lbs = cc.lbs();
	CoinPackedVector ubs = cc.ubs();
	for (int i = 0; i < lbs.getNumElements(); ++i) {
		int col = lbs.getIndices()[i];
		double lb = lbs.getElements()[i];
		if (lb > clbd_[col])
		{
			SCIPchgVarLb(scip_, vars_[col], lb);
			clbd_[col] = lb;
		}
	}
	for (int i = 0; i < ubs.getNumElements(); ++i) {
		int col = ubs.getIndices()[i];
		double ub = ubs.getElements()[i];
		if (ub > cubd_[col]) 
		{
			SCIPchgVarUb(scip_, vars_[col], ub);
			cubd_[col] = ub;
		}
	}
}
