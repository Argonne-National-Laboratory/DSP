/*
 * OsiOoqpSolverInterface.cpp
 *
 *  Created on: Sep 9, 2016
 *      Author: kibaekkim
 */

/** Coin */
#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
/** Dsp */
#include "SolverInterface/OsiOoqpSolverInterface.h"
/** Ooqp */
#include "Status.h"
#include "QpGenSparseMa27.h"
#include "GondzioSolver.h"
#include "MehrotraSolver.h"

// Free memory pointer to by a int pointer
inline void freeCacheInt(int*& ptr) {
	if (ptr != NULL) {
		delete[] ptr;
		ptr = NULL;
	}
}

// Free memory pointer to by a double pointer
inline void freeCacheDouble(double*& ptr) {
	if (ptr != NULL) {
		delete[] ptr;
		ptr = NULL;
	}
}

// Free memory pointer to by a char pointer
inline void freeCacheChar(char*& ptr) {
	if (ptr != NULL) {
		delete[] ptr;
		ptr = NULL;
	}
}

// Free memory pointer to by a CoinPackedMatrix pointer
inline void freeCacheMatrix(CoinPackedMatrix*& ptr) {
	if (ptr != NULL) {
		delete ptr;
		ptr = NULL;
	}
}

void OsiOoqpSolverInterface::freeCachedResults() {
	freeCacheDouble(x_);
	freeCacheDouble(price_);
	freeCacheDouble(reduced_);
	freeCacheDouble(activity_);
}

void OsiOoqpSolverInterface::freeCachedData() {
	freeCacheDouble(dQ_);
	freeCacheMatrix(mat_);
	freeCacheDouble(clbd_);
	freeCacheDouble(cubd_);
	freeCacheDouble(obj_);
	freeCacheDouble(rlbd_);
	freeCacheDouble(rubd_);
	freeCacheChar(sense_);
	freeCacheDouble(rhs_);
	freeCacheDouble(range_);
	nnzQ_ = 0;
	freeCacheInt(irowQ_);
	freeCacheInt(jcolQ_);
}

void OsiOoqpSolverInterface::freeCachedOoqp() {
	if (qpgen_ != NULL) {
		delete qpgen_;
		qpgen_ = NULL;
	}
	if (prob_ != NULL) {
		delete prob_;
		prob_ = NULL;
	}
	if (vars_ != NULL) {
		delete vars_;
		vars_ = NULL;
	}
	if (resid_ != NULL) {
		delete resid_;
		resid_ = NULL;
	}
	if (solver_ != NULL) {
		delete solver_;
		solver_ = NULL;
	}
}

void OsiOoqpSolverInterface::initializeProblem(const CoinPackedMatrix& matrix) {
	/** TODO Quadratic term */
	mat_   = new CoinPackedMatrix(matrix);
	clbd_  = new double [mat_->getNumCols()];
	cubd_  = new double [mat_->getNumCols()];
	obj_   = new double [mat_->getNumCols()];
	rlbd_  = new double [mat_->getNumRows()];
	rubd_  = new double [mat_->getNumRows()];
	sense_ = new char [mat_->getNumRows()];
	rhs_   = new double [mat_->getNumRows()];
	range_ = new double [mat_->getNumRows()];
}

void OsiOoqpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
		const double* collb, const double* colub, const double* obj, int nnzQ,
		const int* irowQ, const int* jcolQ, const double* dQ,
		const double* rowlb, const double* rowub) {
	loadProblem(matrix, collb, colub, obj, rowlb, rowub);
	setQuadraticObjective(nnzQ, irowQ, jcolQ, dQ);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
		double*& collb, double*& colub, double*& obj, int nnzQ, int*& irowQ,
		int*& jcolQ, double*& dQ, double*& rowlb, double*& rowub) {
	loadProblem(*matrix, collb, colub, obj, nnzQ, irowQ, jcolQ, dQ, rowlb, rowub);
	freeCacheMatrix(matrix);
	freeCacheDouble(collb);
	freeCacheDouble(colub);
	freeCacheDouble(obj);
	freeCacheInt(irowQ);
	freeCacheInt(jcolQ);
	freeCacheDouble(dQ);
	freeCacheDouble(rowlb);
	freeCacheDouble(rowub);
}

void OsiOoqpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
		const double* collb, const double* colub, const double* obj, int nnzQ,
		const int* irowQ, const int* jcolQ, const double* dQ,
		const char* rowsen, const double* rowrhs, const double* rowrng) {
	loadProblem(matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
	setQuadraticObjective(nnzQ, irowQ, jcolQ, dQ);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
		double*& collb, double*& colub, double*& obj, int nnzQ, int*& irowQ,
		int*& jcolQ, double*& dQ, char*& rowsen, double*& rowrhs,
		double*& rowrng) {
	loadProblem(*matrix, collb, colub, obj, nnzQ, irowQ, jcolQ, dQ, rowsen, rowrhs, rowrng);
	freeCacheMatrix(matrix);
	freeCacheDouble(collb);
	freeCacheDouble(colub);
	freeCacheDouble(obj);
	freeCacheInt(irowQ);
	freeCacheInt(jcolQ);
	freeCacheDouble(dQ);
	freeCacheChar(rowsen);
	freeCacheDouble(rowrhs);
	freeCacheDouble(rowrng);
}

void OsiOoqpSolverInterface::loadProblem(const int numcols, const int numrows,
		const CoinBigIndex* start, const int* index, const double* value,
		const double* collb, const double* colub, const double* obj, int nnzQ,
		const int* irowQ, const int* jcolQ, const double* dQ,
		const double* rowlb, const double* rowub) {
	/** column length */
	int* len = new int [numcols];
	for (int j = 0; j < numcols; ++j)
		len[j] = start[j+1] - start[j];
	/** create matrix */
	CoinPackedMatrix* matrix = new CoinPackedMatrix(true, numrows, numcols, start[numcols], value, index, start, len);
	/** load problem */
	loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
	/** release memory */
	freeCacheMatrix(matrix);
	freeCacheInt(len);
	/** set quadratic objective function */
	setQuadraticObjective(nnzQ, irowQ, jcolQ, dQ);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::loadProblem(const int numcols, const int numrows,
		const CoinBigIndex* start, const int* index, const double* value,
		const double* collb, const double* colub, const double* obj, int nnzQ,
		const int* irowQ, const int* jcolQ, const double* dQ,
		const char* rowsen, const double* rowrhs, const double* rowrng) {
	/** column length */
	int* len = new int [numcols];
	for (int j = 0; j < numcols; ++j)
		len[j] = start[j+1] - start[j];
	/** create matrix */
	CoinPackedMatrix* matrix = new CoinPackedMatrix(true, numrows, numcols, start[numcols], value, index, start, len);
	/** load problem */
	loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
	/** release memory */
	freeCacheMatrix(matrix);
	freeCacheInt(len);
	/** set quadratic objective function */
	setQuadraticObjective(nnzQ, irowQ, jcolQ, dQ);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setQuadraticObjective(
		int nnz,            /**< [in] number of nonzeros */
		const int* irow,    /**< [in] row indices */
		const int* jcol,    /**< [in] column indices */
		const double* value /**< [in] matrix elements */) {
	nnzQ_ = nnz;
	irowQ_ = new int [nnzQ_];
	jcolQ_ = new int [nnzQ_];
	dQ_ = new double [nnzQ_];
	CoinCopyN(irow, nnz, irowQ_);
	CoinCopyN(jcol, nnz, jcolQ_);
	CoinCopyN(value, nnz, dQ_);
}

void OsiOoqpSolverInterface::writeMps(const char* filename,
		const char* extension, double objSense) const {
	throw CoinError("Not supported.", "writeMps", "OsiOoqpSolverInterface");
}

void OsiOoqpSolverInterface::convertOsiToOoqp(QpGen*& qpgen, QpGenData*& prob) {
	int nx = mat_->getNumCols();
	int my = 0; /**< number of equality constraints */
	int mz = 0; /**< number of inequality constraints */
	int nnzA = 0; /**< number of nonzeros in the equality constraint matrix */
	int nnzC = 0; /**< number of nonzeros in the inequality constraint matrix */

	/** row-wise matrix */
	const CoinPackedMatrix* mat = getMatrixByRow();

	/** count nonzeros */
	for (int i = 0; i < mat->getNumRows(); ++i) {
		const CoinShallowPackedVector row = mat->getVector(i);
		if (sense_[i] == 'E') {
			my++;
			nnzA += row.getNumElements();
		} else {
			mz++;
			nnzC += row.getNumElements();
		}
	}

	/** create QP generator */
	qpgen = new QpGenSparseMa27(nx, my, mz, nnzQ_, nnzA, nnzC);
	//qpgen = new QpGenSparseMa57(nx, my, mz, nnzQ_, nnzA, nnzC);

	/** allocate memory */
	double * xlow  = new double [nx]; /**< variable lower bounds */
	char *   ixlow = new char [nx];   /**< 1 if variable lower bound exists; 0 otherwise. */
	double * xupp  = new double [nx]; /**< variable upper bounds */
	char *   ixupp = new char [nx];   /**< 1 if variable upper bound exists; 0 otherwise. */
	int * irowA    = new int [nnzA];
	int * jcolA    = new int [nnzA];
	double * dA    = new double [nnzA];
	double * b     = new double [my];
	int * irowC    = new int [nnzC];
	int * jcolC    = new int [nnzC];
	double * dC    = new double [nnzC];
	double * clow  = new double [mz];
	char * iclow   = new char [mz];
	double * cupp  = new double [mz];
	char * icupp   = new char [mz];

	double inf = getInfinity();
	for (int j = 0; j < mat->getNumCols(); ++j) {
		xlow[j] = clbd_[j];
		ixlow[j] = 1;
		xupp[j] = cubd_[j];
		ixupp[j] = 1;
		if (clbd_[j] <= -inf) {
			xlow[j] = 0.0;
			ixlow[j] = 0;
		}
		if (cubd_[j] >= inf) {
			xupp[j] = 0.0;
			ixupp[j] = 0;
		}
	}

	int posA = 0;
	int posC = 0;
	for (int i = 0; i < mat->getNumRows(); ++i) {
		const CoinShallowPackedVector row = mat->getVector(i);
		if (sense_[i] == 'E') {
			for (int j = 0; j < row.getNumElements(); ++j) {
				irowA[posA] = i;
				jcolA[posA] = row.getIndices()[j];
				dA[posA] = row.getElements()[j];
				posA++;
			}
			b[i] = rhs_[i];
		} else {
			for (int j = 0; j < row.getNumElements(); ++j) {
				irowC[posC] = i;
				jcolC[posC] = row.getIndices()[j];
				dC[posC] = row.getElements()[j];
				posC++;
			}
			clow[i] = rlbd_[i];
			iclow[i] = 1;
			cupp[i] = rubd_[i];
			icupp[i] = 1;
			if (sense_[i] == 'L') {
				clow[i] = 0.0;
				iclow[i] = 0;
			} else if (sense_[i] == 'G') {
				cupp[i] = 0.0;
				icupp[i] = 0;
			}
		}
	}

	prob = (QpGenData*)dynamic_cast<QpGenSparseSeq*>(qpgen)->copyDataFromSparseTriple(
			obj_,  irowQ_, nnzQ_, jcolQ_, dQ_,
			xlow,  ixlow,  xupp,  ixupp,
			irowA, nnzA,   jcolA, dA,     b,
			irowC, nnzC,   jcolC, dC,
			clow,  iclow,  cupp,  icupp);

	/** release memory */
	freeCacheDouble(xlow);
	freeCacheChar(ixlow);
	freeCacheDouble(xupp);
	freeCacheChar(ixupp);
	freeCacheInt(irowA);
	freeCacheInt(jcolA);
	freeCacheDouble(dA);
	freeCacheDouble(b);
	freeCacheInt(irowC);
	freeCacheInt(jcolC);
	freeCacheDouble(dC);
	freeCacheDouble(clow);
	freeCacheChar(iclow);
	freeCacheDouble(cupp);
	freeCacheChar(icupp);
}

void OsiOoqpSolverInterface::initialSolve() {

	/** free Ooqp objects */
	freeCachedOoqp();

	/** convert Osi data to Ooqp data */
	convertOsiToOoqp(qpgen_, prob_);
	//prob_->print();

	/** declare variables */
	vars_ = (QpGenVars*)qpgen_->makeVariables(prob_);

	/** declare residuals */
	resid_ = (QpGenResiduals*)qpgen_->makeResiduals(prob_);

	/** create solver */
	//solver_ = new GondzioSolver(qp_, prob_);
	solver_ = new MehrotraSolver(qpgen_, prob_);

	/** actual solve */
	gutsOfSolve();
}

void OsiOoqpSolverInterface::resolve() {
	if (updated_)
		initialSolve();
	else
		gutsOfSolve();
}

void OsiOoqpSolverInterface::gutsOfSolve() {
	/** solve */
	status_ = solver_->solve(prob_, vars_, resid_);
	if (handler_->logLevel() >= 4)
		solver_->monitorSelf();

	if (status_ == SUCCESSFUL_TERMINATION || status_ == MAX_ITS_EXCEEDED) {
		/** allocate memory for OOQP solutions */
		double * y      = new double [prob_->my]; /**< dual variables corresponding to equality constraints */
		double * lambda = new double [prob_->mz]; /**< dual variables corresponding to inequality (>=) constraints */
		double * pi     = new double [prob_->mz]; /**< dual variables corresponding to inequality (<=) constraints */
		double * gamma  = new double [prob_->nx]; /**< dual variables corresponding to column lower bounds */
		double * phi    = new double [prob_->nx]; /**< dual variables corresponding to column upper bounds */

		/** allocate memory for Osi solutions */
		freeCachedResults();
		x_ = new double [prob_->nx];
		price_ = new double [mat_->getNumRows()];
		reduced_ = new double [mat_->getNumRows()];
		activity_ = new double [mat_->getNumRows()];

		/** retrieve OOQP solutions */
		objval_ = prob_->objectiveValue(vars_);
		vars_->x->copyIntoArray(x_);
		vars_->y->copyIntoArray(y);
		vars_->lambda->copyIntoArray(lambda);
		vars_->pi->copyIntoArray(pi);
		vars_->gamma->copyIntoArray(gamma);
		vars_->phi->copyIntoArray(phi);
		nIters_ = solver_->iter;

		/** convert OOQP solutions to Osi solutions */
		int pos1 = 0, pos2 = 0, pos3 = 0;
		for (int i = 0; i < mat_->getNumRows(); ++i) {
			if (sense_[i] == 'E')
				price_[i] = y[pos1++];
			else if (sense_[i] == 'G')
				price_[i] = lambda[pos2++];
			else if (sense_[i] == 'L')
				price_[i] = pi[pos3++];
			else if (sense_[i] == 'R')
				price_[i] = lambda[pos2++] - pi[pos3++];
		}
		for (int j = 0; j < mat_->getNumCols(); ++j) {
			reduced_[j] = gamma[j] - phi[j];
		}
		/** TODO: Implement assigning activity */

		/** release memory */
		freeCacheDouble(y);
		freeCacheDouble(lambda);
		freeCacheDouble(pi);
		freeCacheDouble(gamma);
		freeCacheDouble(phi);
	}
	/** say the model is old. */
	updated_ = false;
}

void OsiOoqpSolverInterface::branchAndBound() {
	throw CoinError("Not supported.", "branchAndBound", "OsiOoqpSolverInterface");
}

bool OsiOoqpSolverInterface::isAbandoned() const {
	return (status_ == NOT_FINISHED);
}

bool OsiOoqpSolverInterface::isProvenOptimal() const {
	return (status_ == SUCCESSFUL_TERMINATION);
}

bool OsiOoqpSolverInterface::isProvenPrimalInfeasible() const {
	/** TODO: How to check primal infeasibility? */
	CoinError("Probably infeasible.", "isProvenPrimalInfeasible", "OsiOoqpSolverInterface");
	return (status_ == INFEASIBLE);
}

bool OsiOoqpSolverInterface::isProvenDualInfeasible() const {
	/** TODO: How to check dual infeasibility? */
	CoinError("Probably dual infeasible.", "isProvenDualInfeasible", "OsiOoqpSolverInterface");
	return (status_ == INFEASIBLE);
}

bool OsiOoqpSolverInterface::isIterationLimitReached() const {
	return (status_ == MAX_ITS_EXCEEDED);
}

CoinWarmStart* OsiOoqpSolverInterface::getEmptyWarmStart() const {
	CoinError("Not supported.", "getEmptyWarmStart", "OsiOoqpSolverInterface");
	return NULL;
}

CoinWarmStart* OsiOoqpSolverInterface::getWarmStart() const {
	CoinError("Not supported.", "getWarmStart", "OsiOoqpSolverInterface");
	return NULL;
}

bool OsiOoqpSolverInterface::setWarmStart(const CoinWarmStart* warmstart) {
	CoinError("Not supported.", "setWarmStart", "OsiOoqpSolverInterface");
	return false;
}

int OsiOoqpSolverInterface::getNumCols() const {
	return mat_->getNumCols();
}

int OsiOoqpSolverInterface::getNumRows() const {
	return mat_->getNumRows();
}

int OsiOoqpSolverInterface::getNumElements() const {
	return mat_->getNumElements();
}

const double* OsiOoqpSolverInterface::getColLower() const {
	return clbd_;
}

const double* OsiOoqpSolverInterface::getColUpper() const {
	return cubd_;
}

const char* OsiOoqpSolverInterface::getRowSense() const {
	return sense_;
}

const double* OsiOoqpSolverInterface::getRightHandSide() const {
	return rhs_;
}

const double* OsiOoqpSolverInterface::getRowRange() const {
	return range_;
}

const double* OsiOoqpSolverInterface::getRowLower() const {
	return rlbd_;
}

const double* OsiOoqpSolverInterface::getRowUpper() const {
	return rubd_;
}

/** TODO: Returns linear objective coefficients only */
const double* OsiOoqpSolverInterface::getObjCoefficients() const {
	return obj_;
}

double OsiOoqpSolverInterface::getObjSense() const {
	return obj_sense_;
}

bool OsiOoqpSolverInterface::isContinuous(int colIndex) const {
	return true;
}

const CoinPackedMatrix* OsiOoqpSolverInterface::getMatrixByRow() const {
	if (mat_->isColOrdered())
		mat_->reverseOrdering();
	return mat_;
}

const CoinPackedMatrix* OsiOoqpSolverInterface::getMatrixByCol() const {
	if (mat_->isColOrdered() == false)
		mat_->reverseOrdering();
	return mat_;
}

double OsiOoqpSolverInterface::getInfinity() const {
	return 1.0e+20;
}

const double* OsiOoqpSolverInterface::getColSolution() const {
	return x_;
}

const double* OsiOoqpSolverInterface::getRowPrice() const {
	return price_;
}

const double* OsiOoqpSolverInterface::getReducedCost() const {
	return reduced_;
}

const double* OsiOoqpSolverInterface::getRowActivity() const {
	return activity_;
}

double OsiOoqpSolverInterface::getObjValue() const {
	return objval_;
}

int OsiOoqpSolverInterface::getIterationCount() const {
	return nIters_;
}

std::vector<double*> OsiOoqpSolverInterface::getDualRays(
		int maxNumRays,
		bool fullRay) const {
	/** TODO: Implement */
	CoinError("Not implemented.", "getDualRays", "OsiOoqpSolverInterface");
	return std::vector<double*>();
}

std::vector<double*> OsiOoqpSolverInterface::getPrimalRays(int maxNumRays) const {
	/** TODO: Implement */
	CoinError("Not implemented.", "getPrimalRays", "OsiOoqpSolverInterface");
	return std::vector<double*>();
}

void OsiOoqpSolverInterface::setObjCoeff(
		int elementIndex,
		double elementValue) {
	if (elementIndex < mat_->getNumCols())
		obj_[elementIndex] = elementValue;
	else
		CoinError("Column index is out of range.", "setObjCoeff", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setObjSense(double s) {
	obj_sense_ = s;
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setColLower(int elementIndex,
		double elementValue) {
	if (elementIndex < mat_->getNumCols())
		clbd_[elementIndex] = elementValue;
	else
		CoinError("Column index is out of range.", "setColLower", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setColUpper(int elementIndex,
		double elementValue) {
	if (elementIndex < mat_->getNumCols())
		cubd_[elementIndex] = elementValue;
	else
		CoinError("Column index is out of range.", "setColUpper", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setRowLower(int elementIndex,
		double elementValue) {
	if (elementIndex < mat_->getNumRows())
		rlbd_[elementIndex] = elementValue;
	else
		CoinError("Row index is out of range.", "setRowLower", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setRowUpper(int elementIndex,
		double elementValue) {
	if (elementIndex < mat_->getNumRows())
		rlbd_[elementIndex] = elementValue;
	else
		CoinError("Row index is out of range.", "setRowUpper", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setRowType(int index, char sense,
		double rightHandSide, double range) {
	double rowlb, rowub;
	convertSenseToBound(sense, rightHandSide, range, rowlb, rowub);
	sense_[index] = sense;
	range_[index] = range;
	rhs_[index] = rightHandSide;
	rlbd_[index] = rowlb;
	rubd_[index] = rowub;
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setColSolution(const double* colsol) {
	/** TODO: Not supported? */
	CoinError("Not supported.", "setColSolution", "OsiOoqpSolverInterface");
}

void OsiOoqpSolverInterface::setRowPrice(const double* rowprice) {
	/** TODO: Not supported */
	CoinError("Not supported.", "setRowPrice", "OsiOoqpSolverInterface");
}

void OsiOoqpSolverInterface::setContinuous(int index) {
	CoinError("Not supported.", "setContinuous", "OsiOoqpSolverInterface");
}

void OsiOoqpSolverInterface::setInteger(int index) {
	CoinError("Not supported.", "setInteger", "OsiOoqpSolverInterface");
}

void OsiOoqpSolverInterface::addCol(
		const CoinPackedVectorBase& vec,
		const double collb,
		const double colub,
		const double obj) {
	/** add a column */
	mat_->appendCol(vec);
	/** store original data */
	double* clbd0 = clbd_;
	double* cubd0 = cubd_;
	double* obj0 = obj_;
	/** reallocate memory */
	clbd_ = NULL;
	cubd_ = NULL;
	obj_ = NULL;
	delete [] x_; x_ = NULL;
	clbd_ = new double [mat_->getNumCols()];
	cubd_ = new double [mat_->getNumCols()];
	obj_  = new double [mat_->getNumCols()];
	x_    = new double [mat_->getNumCols()];
	/** copy original data */
	CoinCopyN(clbd0, mat_->getNumCols()-1, clbd_);
	CoinCopyN(cubd0, mat_->getNumCols()-1, cubd_);
	CoinCopyN(obj0, mat_->getNumCols()-1, obj_);
	/** add new data */
	clbd_[mat_->getNumCols()-1] = collb;
	cubd_[mat_->getNumCols()-1] = colub;
	obj_[mat_->getNumCols()-1] = obj;
	/** release original memory */
	freeCacheDouble(clbd0);
	freeCacheDouble(cubd0);
	freeCacheDouble(obj0);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::deleteCols(const int num, const int* colIndices) {
	/** number of columns before deleted */
	int nints0 = mat_->getNumCols();
	/** delete columns */
	mat_->deleteCols(num, colIndices);
	/** store original data */
	double* clbd0 = clbd_;
	double* cubd0 = cubd_;
	double* obj0 = obj_;
	/** reallocate memory */
	clbd_ = NULL;
	cubd_ = NULL;
	obj_ = NULL;
	delete [] x_; x_ = NULL;
	clbd_ = new double [mat_->getNumCols()];
	cubd_ = new double [mat_->getNumCols()];
	obj_  = new double [mat_->getNumCols()];
	x_    = new double [mat_->getNumCols()];

	/** column indices to delete */
	std::vector<int> delCols;
	for (int j = 0; j < num; ++j)
		delCols.push_back(colIndices[j]);
	/** sort the indices */
	std::sort(delCols.begin(), delCols.end());

	/** copy data */
	for (int j = 0, k = 0; j < nints0; ++j) {
		if (delCols[k] == j) {
			k++;
			continue;
		}
		clbd_[j-k] = clbd0[j];
		cubd_[j-k] = cubd0[j];
		obj_[j-k] = obj0[j];
	}

	/** release original memory */
	freeCacheDouble(clbd0);
	freeCacheDouble(cubd0);
	freeCacheDouble(obj0);

	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::addRow(
		const CoinPackedVectorBase& vec,
		const double rowlb,
		const double rowub) {
	/** append row */
	mat_->appendRow(vec);
	/** store original data */
	double* rlbd0 = rlbd_;
	double* rubd0 = rubd_;
	char* sense0 = sense_;
	double* rhs0 = rhs_;
	double* range0 = range_;
	/** reallocate memory */
	rlbd_ = NULL;
	rubd_ = NULL;
	sense_ = NULL;
	rhs_ = NULL;
	range_ = NULL;
	delete [] price_; price_ = NULL;
	delete [] reduced_; reduced_ = NULL;
	delete [] activity_; activity_ = NULL;
	rlbd_ = new double [mat_->getNumRows()];
	rubd_ = new double [mat_->getNumRows()];
	sense_ = new char [mat_->getNumRows()];
	rhs_ = new double [mat_->getNumRows()];
	range_ = new double [mat_->getNumRows()];
	price_ = new double [mat_->getNumRows()];
	reduced_  = new double [mat_->getNumRows()];
	activity_ = new double [mat_->getNumRows()];
	/** copy original data */
	CoinCopyN(rlbd0, mat_->getNumRows()-1, rlbd_);
	CoinCopyN(rubd0, mat_->getNumRows()-1, rubd_);
	CoinCopyN(sense0, mat_->getNumRows()-1, sense_);
	CoinCopyN(rhs0, mat_->getNumRows()-1, rhs_);
	CoinCopyN(range0, mat_->getNumRows()-1, range_);
	/** add new data */
	rlbd_[mat_->getNumRows()-1] = rowlb;
	rubd_[mat_->getNumRows()-1] = rowub;
	convertBoundToSense(rowlb, rowub, sense_[mat_->getNumRows()-1],
			rhs_[mat_->getNumRows()-1], range_[mat_->getNumRows()-1]);
	/** release original memory */
	freeCacheDouble(rlbd0);
	freeCacheDouble(rubd0);
	freeCacheChar(sense0);
	freeCacheDouble(rhs0);
	freeCacheDouble(range0);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::addRow(
		const CoinPackedVectorBase& vec,
		const char rowsen,
		const double rowrhs,
		const double rowrng) {
	double rowlb, rowub;
	convertSenseToBound(rowsen, rowrhs, rowrng, rowlb, rowub);
	addRow(vec, rowlb, rowub);
}

void OsiOoqpSolverInterface::deleteRows(const int num, const int* rowIndices) {
	/** number of rows before deleted */
	int nrows0 = mat_->getNumRows();

	/** delete rows */
	mat_->deleteRows(num, rowIndices);

	/** store original data */
	double* rlbd0 = rlbd_;
	double* rubd0 = rubd_;
	char* sense0 = sense_;
	double* rhs0 = rhs_;
	double* range0 = range_;

	/** reallocate memory */
	rlbd_ = NULL;
	rubd_ = NULL;
	sense_ = NULL;
	rhs_ = NULL;
	range_ = NULL;
	delete [] price_; price_ = NULL;
	delete [] reduced_; reduced_ = NULL;
	delete [] activity_; activity_ = NULL;
	rlbd_ = new double [mat_->getNumRows()];
	rubd_ = new double [mat_->getNumRows()];
	sense_ = new char [mat_->getNumRows()];
	rhs_ = new double [mat_->getNumRows()];
	range_ = new double [mat_->getNumRows()];
	price_ = new double [mat_->getNumRows()];
	reduced_  = new double [mat_->getNumRows()];
	activity_ = new double [mat_->getNumRows()];

	/** row indices to delete */
	std::vector<int> delRows;
	for (int j = 0; j < num; ++j)
		delRows.push_back(rowIndices[j]);
	/** sort the indices */
	std::sort(delRows.begin(), delRows.end());

	/** copy data */
	for (int j = 0, k = 0; j < nrows0; ++j) {
		if (delRows[k] == j) {
			k++;
			continue;
		}
		rlbd_[j-k] = rlbd0[j];
		rubd_[j-k] = rubd0[j];
		sense_[j-k] = sense0[j];
		rhs_[j-k] = rhs0[j];
		range_[j-k] = range0[j];
	}

	/** release original memory */
	freeCacheDouble(rlbd0);
	freeCacheDouble(rubd0);
	freeCacheChar(sense0);
	freeCacheDouble(rhs0);
	freeCacheDouble(range0);

	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
		const double* collb, const double* colub, const double* obj,
		const double* rowlb, const double* rowub) {
	freeCachedResults();
	freeCachedData();
	/** allocate memory */
	initializeProblem(matrix);
	/** copy data */
	CoinCopyN(collb, mat_->getNumCols(), clbd_);
	CoinCopyN(colub, mat_->getNumCols(), cubd_);
	CoinCopyN(obj,   mat_->getNumCols(), obj_);
	CoinCopyN(rowlb, mat_->getNumRows(), rlbd_);
	CoinCopyN(rowub, mat_->getNumRows(), rubd_);
	for (int i = 0; i < mat_->getNumRows(); ++i)
		convertBoundToSense(rowlb[i], rowub[i], sense_[i], rhs_[i], range_[i]);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
		double*& collb, double*& colub, double*& obj, double*& rowlb,
		double*& rowub) {
	loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
	freeCacheMatrix(matrix);
	freeCacheDouble(collb);
	freeCacheDouble(colub);
	freeCacheDouble(obj);
	freeCacheDouble(rowlb);
	freeCacheDouble(rowub);
}

void OsiOoqpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
		const double* collb, const double* colub, const double* obj,
		const char* rowsen, const double* rowrhs, const double* rowrng) {
	freeCachedResults();
	freeCachedData();
	/** allocate memory */
	initializeProblem(matrix);
	/** copy data */
	CoinCopyN(collb, mat_->getNumCols(), clbd_);
	CoinCopyN(colub, mat_->getNumCols(), cubd_);
	CoinCopyN(obj,   mat_->getNumCols(), obj_);
	CoinCopyN(rowsen, mat_->getNumRows(), sense_);
	CoinCopyN(rowrhs, mat_->getNumRows(), rhs_);
	CoinCopyN(rowrng, mat_->getNumRows(), range_);
	for (int i = 0; i < mat_->getNumRows(); ++i)
		convertSenseToBound(sense_[i], rhs_[i], range_[i], rlbd_[i], rubd_[i]);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
		double*& collb, double*& colub, double*& obj, char*& rowsen,
		double*& rowrhs, double*& rowrng) {
	loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
	freeCacheMatrix(matrix);
	freeCacheDouble(collb);
	freeCacheDouble(colub);
	freeCacheDouble(obj);
	freeCacheChar(rowsen);
	freeCacheDouble(rowrhs);
	freeCacheDouble(rowrng);
}

void OsiOoqpSolverInterface::loadProblem(const int numcols, const int numrows,
		const CoinBigIndex* start, const int* index, const double* value,
		const double* collb, const double* colub, const double* obj,
		const double* rowlb, const double* rowub) {
	/** column length */
	int* len = new int [numcols];
	for (int j = 0; j < numcols; ++j)
		len[j] = start[j+1] - start[j];
	/** create matrix */
	CoinPackedMatrix* matrix = new CoinPackedMatrix(true, numrows, numcols, start[numcols], value, index, start, len);
	/** load problem */
	loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
	/** release memory */
	freeCacheMatrix(matrix);
	freeCacheInt(len);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::loadProblem(const int numcols, const int numrows,
		const CoinBigIndex* start, const int* index, const double* value,
		const double* collb, const double* colub, const double* obj,
		const char* rowsen, const double* rowrhs, const double* rowrng) {
	/** column length */
	int* len = new int [numcols];
	for (int j = 0; j < numcols; ++j)
		len[j] = start[j+1] - start[j];
	/** create matrix */
	CoinPackedMatrix* matrix = new CoinPackedMatrix(true, numrows, numcols, start[numcols], value, index, start, len);
	/** load problem */
	loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
	/** release memory */
	freeCacheMatrix(matrix);
	freeCacheInt(len);
	/** say the model is updated. */
	updated_ = true;
}

OsiOoqpSolverInterface::OsiOoqpSolverInterface() :
	obj_sense_(1.0),
	mat_(NULL),
	clbd_(NULL),
	cubd_(NULL),
	obj_(NULL),
	rlbd_(NULL),
	rubd_(NULL),
	sense_(NULL),
	rhs_(NULL),
	range_(NULL),
	status_(UNKNOWN),
	objval_(0.0),
	x_(NULL),
	price_(NULL),
	reduced_(NULL),
	activity_(NULL),
	nIters_(0),
	qpgen_(NULL),
	prob_(NULL),
	vars_(NULL),
	resid_(NULL),
	solver_(NULL),
	nnzQ_(0),
	irowQ_(NULL),
	jcolQ_(NULL),
	dQ_(NULL),
	updated_(true) {
	/**< nothing to do */
}

OsiOoqpSolverInterface* OsiOoqpSolverInterface::clone(bool copyData) const {
	return (new OsiOoqpSolverInterface(*this));
}

OsiOoqpSolverInterface::OsiOoqpSolverInterface(const OsiOoqpSolverInterface& rhs) :
	OsiSolverInterface(rhs),
	obj_sense_(rhs.obj_sense_),
	status_(rhs.status_),
	objval_(rhs.objval_),
	nIters_(rhs.nIters_),
	qpgen_(NULL),
	prob_(NULL),
	vars_(NULL),
	resid_(NULL),
	solver_(NULL),
	updated_(true) {
	/** load problem */
	loadProblem(*(rhs.mat_), rhs.clbd_, rhs.cubd_, rhs.obj_,
			rhs.nnzQ_, rhs.irowQ_, rhs.jcolQ_, rhs.dQ_, rhs.rlbd_, rhs.rubd_);
	/** copy solution data */
	x_ = new double [mat_->getNumCols()];
	price_ = new double [mat_->getNumRows()];
	reduced_ = new double [mat_->getNumRows()];
	activity_ = new double [mat_->getNumRows()];
	CoinCopyN(rhs.x_, mat_->getNumCols(), x_);
	CoinCopyN(rhs.price_, mat_->getNumRows(), price_);
	CoinCopyN(rhs.reduced_, mat_->getNumRows(), reduced_);
	CoinCopyN(rhs.activity_, mat_->getNumRows(), activity_);
}

OsiOoqpSolverInterface& OsiOoqpSolverInterface::operator =(
		const OsiOoqpSolverInterface& rhs) {
	if (this != &rhs) {
		/** free Ooqp */
		freeCachedOoqp();
		/** copy */
		obj_sense_ = rhs.obj_sense_;
		status_ = rhs.status_;
		objval_ = rhs.objval_;
		nIters_ = rhs.nIters_;
		updated_ = rhs.updated_;
		/** load problem */
		loadProblem(*(rhs.mat_), rhs.clbd_, rhs.cubd_, rhs.obj_,
				rhs.nnzQ_, rhs.irowQ_, rhs.jcolQ_, rhs.dQ_, rhs.rlbd_, rhs.rubd_);
		/** copy solution data */
		x_ = new double [mat_->getNumCols()];
		price_ = new double [mat_->getNumRows()];
		reduced_ = new double [mat_->getNumRows()];
		activity_ = new double [mat_->getNumRows()];
		CoinCopyN(rhs.x_, mat_->getNumCols(), x_);
		CoinCopyN(rhs.price_, mat_->getNumRows(), price_);
		CoinCopyN(rhs.reduced_, mat_->getNumRows(), reduced_);
		CoinCopyN(rhs.activity_, mat_->getNumRows(), activity_);
	}
	return *this;
}

OsiOoqpSolverInterface::~OsiOoqpSolverInterface() {
	freeCachedResults();
	freeCachedOoqp();
	freeCachedData();
}

void OsiOoqpSolverInterface::applyRowCut(const OsiRowCut& rc) {
	addRow(rc.row(), rc.lb(), rc.ub());
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::applyColCut(const OsiColCut& cc) {
	CoinPackedVector lbs = cc.lbs();
	CoinPackedVector ubs = cc.ubs();
	for (int i = 0; i < lbs.getNumElements(); ++i) {
		int col = lbs.getIndices()[i];
		double lb = lbs.getElements()[i];
		if (lb > clbd_[col]) clbd_[col] = lb;
	}
	for (int i = 0; i < ubs.getNumElements(); ++i) {
		int col = ubs.getIndices()[i];
		double ub = ubs.getElements()[i];
		if (ub > cubd_[col]) cubd_[col] = ub;
	}
	/** say the model is updated. */
	updated_ = true;
}
