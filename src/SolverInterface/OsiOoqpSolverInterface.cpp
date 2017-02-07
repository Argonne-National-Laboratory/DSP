/*
 * OsiOoqpSolverInterface.cpp
 *
 *  Created on: Sep 9, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** Coin */
#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
/** Dsp */
#include "DspMessage.h"
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
	x_.clear();
	price_.clear();
	reduced_.clear();
	activity_.clear();
}

void OsiOoqpSolverInterface::freeCachedData() {
	freeCacheMatrix(mat_);
	clbd_.clear();
	cubd_.clear();
	obj_.clear();
	rlbd_.clear();
	rubd_.clear();
	sense_.clear();
	rhs_.clear();
	range_.clear();
	nnzQ_ = 0;
	irowQ_.clear();
	jcolQ_.clear();
	dQ_.clear();
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
	mat_ = new CoinPackedMatrix(matrix);
	clbd_.resize(mat_->getNumCols());
	cubd_.resize(mat_->getNumCols());
	obj_.resize(mat_->getNumCols());
	rlbd_.resize(mat_->getNumRows());
	rubd_.resize(mat_->getNumRows());
	sense_.resize(mat_->getNumRows());
	rhs_.resize(mat_->getNumRows());
	range_.resize(mat_->getNumRows());
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
	irowQ_.assign(irow, irow + nnz);
	jcolQ_.assign(jcol, jcol + nnz);
	dQ_.assign(value, value + nnz);
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
	DSPdebug(mat->verifyMtx(4));

	/** count nonzeros */
	my_empty_.clear();
	mz_empty_.clear();
	for (int i = 0; i < mat->getNumRows(); ++i) {
		const CoinShallowPackedVector row = mat->getVector(i);
		if (sense_[i] == 'E') {
			if (row.getNumElements() > 0) {
				my++;
				nnzA += row.getNumElements();
			} else
				my_empty_.push_back(i);
		} else {
			if (row.getNumElements() > 0) {
				mz++;
				nnzC += row.getNumElements();
			} else
				mz_empty_.push_back(i);
		}
	}
#ifdef DSP_DEBUG_MORE
	printf("nx(# of cols) %d, my(# of Eq.) %d, mz (# of Ineq.) %d, nnzQ_ %d, nnzA %d, nnzC %d, empty rows %u\n",
			nx, my, mz, nnzQ_, nnzA, nnzC, my_empty_.size() + mz_empty_.size());
#endif

	/** create QP generator */
	qpgen = new QpGenSparseMa27(nx, my, mz, nnzQ_, nnzA, nnzC);
	//qpgen = new QpGenSparseMa57(nx, my, mz, nnzQ_, nnzA, nnzC);

	std::vector<double> xlow(nx, 0.0);
	std::vector<char> ixlow(nx, 0);
	std::vector<double> xupp(nx, 0.0);
	std::vector<char> ixupp(nx, 0);
	std::vector<int> irowA;
	std::vector<int> jcolA;
	std::vector<double> dA;
	std::vector<double> b;
	std::vector<int> irowC;
	std::vector<int> jcolC;
	std::vector<double> dC;
	std::vector<double> clow;
	std::vector<char> iclow;
	std::vector<double> cupp;
	std::vector<char> icupp;

	/** reserve memory */
	irowA.reserve(nnzA);
	jcolA.reserve(nnzA);
	dA.reserve(nnzA);
	b.reserve(my);
	irowC.reserve(nnzC);
	jcolC.reserve(nnzC);
	dC.reserve(nnzC);
	clow.reserve(mz);
	iclow.reserve(mz);
	cupp.reserve(mz);
	icupp.reserve(mz);

	double inf = getInfinity();
	for (int j = 0; j < nx; ++j) {
		if (clbd_[j] > -inf) {
			xlow[j] = clbd_[j];
			ixlow[j] = 1;
		}
		if (cubd_[j] < inf) {
			xupp[j] = cubd_[j];
			ixupp[j] = 1;
		}
#ifdef DSP_DEBUG_MORE
		printf("j %d, xlow %e, ixlow %u, clbd_ %e, xupp %e, ixupp %u, cubd_ %e\n",
				j, xlow[j], ixlow[j], clbd_[j], xupp[j], ixupp[j], cubd_[j]);
#endif
	}

	int iposA = 0, iposC = 0;
	for (int i = 0; i < mat->getNumRows(); ++i) {
		const CoinShallowPackedVector row = mat->getVector(i);
		if (row.getNumElements() == 0) continue;
		if (sense_[i] == 'E') {
			for (int j = 0; j < row.getNumElements(); ++j) {
				irowA.push_back(iposA);
				jcolA.push_back(row.getIndices()[j]);
				dA.push_back(row.getElements()[j]);
#ifdef DSP_DEBUG_MORE
				printf("irowA %d, jcolA %d, dA %e\n", iposA, row.getIndices()[j], row.getElements()[j]);
#endif
			}
			b.push_back(rhs_[i]);
			iposA++;
#ifdef DSP_DEBUG_MORE
			printf("b %e\n", rhs_[i]);
#endif
		} else {
			for (int j = 0; j < row.getNumElements(); ++j) {
				irowC.push_back(iposC);
				jcolC.push_back(row.getIndices()[j]);
				dC.push_back(row.getElements()[j]);
			}
			if (sense_[i] == 'L') {
				clow.push_back(0.0);
				iclow.push_back(0);
				cupp.push_back(rubd_[i]);
				icupp.push_back(1);
			} else if (sense_[i] == 'G') {
				clow.push_back(rlbd_[i]);
				iclow.push_back(1);
				cupp.push_back(0.0);
				icupp.push_back(0);
			} else if (sense_[i] == 'R') {
				clow.push_back(rlbd_[i]);
				iclow.push_back(1);
				cupp.push_back(rubd_[i]);
				icupp.push_back(1);
			}
			iposC++;
#ifdef DSP_DEBUG_MORE
			int endpos = clow.size()-1;
			printf("i %d, clow %e, iclow %u, rlbd_ %e, cupp %e, icupp %u, rubd_ %e\n",
					i, clow[endpos], iclow[endpos], rlbd_[i], cupp[endpos], icupp[endpos], rubd_[i]);
#endif
		}
	}
#ifdef DSP_DEBUG_MORE
	printf("irowA.size() %u irowC.size() %u clow.size() %u\n", irowA.size(), irowC.size(), clow.size());
	DspMessage::printArray(nx, &obj_[0]);
#endif

	prob = (QpGenData*)dynamic_cast<QpGenSparseSeq*>(qpgen)->copyDataFromSparseTriple(
			&obj_[0],  &irowQ_[0], nnzQ_,     &jcolQ_[0], &dQ_[0],
			&xlow[0],  &ixlow[0],  &xupp[0],  &ixupp[0],
			&irowA[0], nnzA,       &jcolA[0], &dA[0],     &b[0],
			&irowC[0], nnzC,       &jcolC[0], &dC[0],
			&clow[0],  &iclow[0],  &cupp[0],  &icupp[0]);
//	prob->print();
}

void OsiOoqpSolverInterface::initialSolve() {

	/** free Ooqp objects */
	freeCachedOoqp();

	/** convert Osi data to Ooqp data */
	convertOsiToOoqp(qpgen_, prob_);

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
	if (messageHandler()->logLevel() >= 4)
		solver_->monitorSelf();
	status_ = solver_->solve(prob_, vars_, resid_);
	DSPdebugMessage("logLevel %d, OOQP status %d\n", messageHandler()->logLevel(), status_);

	if (status_ == SUCCESSFUL_TERMINATION || status_ == MAX_ITS_EXCEEDED) {
		/** allocate memory for OOQP solutions */
		std::vector<double> y(prob_->my);      /**< dual variables corresponding to equality constraints */
		std::vector<double> lambda(prob_->mz); /**< dual variables corresponding to inequality (>=) constraints */
		std::vector<double> pi(prob_->mz);     /**< dual variables corresponding to inequality (<=) constraints */
		std::vector<double> gamma(prob_->nx);  /**< dual variables corresponding to column lower bounds */
		std::vector<double> phi(prob_->nx);    /**< dual variables corresponding to column upper bounds */

		/** allocate memory for Osi solutions */
		freeCachedResults();
		x_.resize(prob_->nx);
		reduced_.resize(mat_->getNumCols());
		activity_.resize(mat_->getNumRows());
		price_.resize(mat_->getNumRows());

		/** retrieve OOQP solutions */
		objval_ = prob_->objectiveValue(vars_);
		vars_->x->copyIntoArray(&x_[0]);
		vars_->y->copyIntoArray(&y[0]);
		vars_->lambda->copyIntoArray(&lambda[0]);
		vars_->pi->copyIntoArray(&pi[0]);
		vars_->gamma->copyIntoArray(&gamma[0]);
		vars_->phi->copyIntoArray(&phi[0]);
		nIters_ = solver_->iter;

		/** convert OOQP solutions to Osi solutions */
		int pos1 = 0, pos2 = 0;
		for (int i = 0, iy = 0, iz = 0; i < mat_->getNumRows(); ++i) {
			price_[i] = 0.0;
			if (my_empty_.size() > iy && my_empty_[iy] == i) {
				iy++;
				continue;
			}
			if (mz_empty_.size() > iz && mz_empty_[iz] == i) {
				iz++;
				continue;
			}
			if (sense_[i] == 'E')
				price_[i] = y[pos1++];
			else if (sense_[i] == 'G')
				price_[i] = lambda[pos2++];
			else if (sense_[i] == 'L')
				price_[i] = -pi[pos2++];
			else if (sense_[i] == 'R')
				price_[i] = lambda[pos2++] - pi[pos2++];
#ifdef DSP_DEBUG_MORE
			printf("sense %c price %e, ", sense_[i], price_[i]);
			switch (sense_[i]) {
			case 'E':
				printf("y %d %e\n", pos1-1, y[pos1-1]);
				break;
			default:
				printf("lambda %d %e pi %e\n", pos2-1, lambda[pos2-1], pi[pos2-1]);
				break;
			}
#endif
		}
#ifdef DSP_DEBUG_MORE
		double pix = 0.0;
		for (int i = 0; i < mat_->getNumRows(); ++i)
			pix += price_[i] * rhs_[i];
		printf("c^T x = %e, pi^T b = %e.\n", objval_, pix);
#endif
		for (int j = 0; j < mat_->getNumCols(); ++j) {
			reduced_[j] = gamma[j] - phi[j];
		}
		/** TODO: Implement assigning activity */
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
	return &clbd_[0];
}

const double* OsiOoqpSolverInterface::getColUpper() const {
	return &cubd_[0];
}

const char* OsiOoqpSolverInterface::getRowSense() const {
	return &sense_[0];
}

const double* OsiOoqpSolverInterface::getRightHandSide() const {
	return &rhs_[0];
}

const double* OsiOoqpSolverInterface::getRowRange() const {
	return &range_[0];
}

const double* OsiOoqpSolverInterface::getRowLower() const {
	return &rlbd_[0];
}

const double* OsiOoqpSolverInterface::getRowUpper() const {
	return &rubd_[0];
}

/** TODO: Returns linear objective coefficients only */
const double* OsiOoqpSolverInterface::getObjCoefficients() const {
	return &obj_[0];
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
	return &x_[0];
}

const double* OsiOoqpSolverInterface::getRowPrice() const {
	return &price_[0];
}

const double* OsiOoqpSolverInterface::getReducedCost() const {
	return &reduced_[0];
}

const double* OsiOoqpSolverInterface::getRowActivity() const {
	return &activity_[0];
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
	if (elementIndex < obj_.size())
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
	if (elementIndex < clbd_.size())
		clbd_[elementIndex] = elementValue;
	else
		CoinError("Column index is out of range.", "setColLower", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setColUpper(int elementIndex,
		double elementValue) {
	if (elementIndex < cubd_.size())
		cubd_[elementIndex] = elementValue;
	else
		CoinError("Column index is out of range.", "setColUpper", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setRowLower(int elementIndex,
		double elementValue) {
	if (elementIndex < rlbd_.size())
		rlbd_[elementIndex] = elementValue;
	else
		CoinError("Row index is out of range.", "setRowLower", "OsiOoqpSolverInterface");
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::setRowUpper(int elementIndex,
		double elementValue) {
	if (elementIndex < rubd_.size())
		rubd_[elementIndex] = elementValue;
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
	clbd_.push_back(collb);
	cubd_.push_back(colub);
	obj_.push_back(obj);
	x_.push_back(0);
	reduced_.push_back(0);
	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::deleteCols(const int num, const int* colIndices) {
	/** number of columns before deleted */
	int nints0 = mat_->getNumCols();

	/** delete columns */
	mat_->deleteCols(num, colIndices);

	/** store original data */
	std::vector<double> clbd0(clbd_);
	std::vector<double> cubd0(cubd_);
	std::vector<double> obj0(obj_);

	/** reallocate memory */
	clbd_.clear();
	cubd_.clear();
	obj_.clear();
	clbd_.reserve(mat_->getNumCols());
	cubd_.reserve(mat_->getNumCols());
	obj_.reserve(mat_->getNumCols());
	x_.resize(mat_->getNumCols());
	reduced_.resize(mat_->getNumCols());

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
		clbd_.push_back(clbd0[j]);
		cubd_.push_back(cubd0[j]);
		obj_.push_back(obj0[j]);
	}

	/** say the model is updated. */
	updated_ = true;
}

void OsiOoqpSolverInterface::addRow(
		const CoinPackedVectorBase& vec,
		const double rowlb,
		const double rowub) {
	/** append row */
	mat_->appendRow(vec);
	rlbd_.push_back(rowlb);
	rubd_.push_back(rowub);
	price_.push_back(0);
	activity_.push_back(0);
	sense_.push_back('L');
	rhs_.push_back(0);
	range_.push_back(0);
	convertBoundToSense(rowlb, rowub, sense_[mat_->getNumRows()-1],
			rhs_[mat_->getNumRows()-1], range_[mat_->getNumRows()-1]);

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
	std::vector<double> rlbd0(rlbd_);
	std::vector<double> rubd0(rubd_);
	std::vector<char> sense0(sense_);
	std::vector<double> rhs0(rhs_);
	std::vector<double> range0(range_);

	/** reallocate memory */
	rlbd_.clear();
	rubd_.clear();
	sense_.clear();
	rhs_.clear();
	range_.clear();
	rlbd_.reserve(mat_->getNumRows());
	rubd_.reserve(mat_->getNumRows());
	sense_.reserve(mat_->getNumRows());
	rhs_.reserve(mat_->getNumRows());
	range_.reserve(mat_->getNumRows());
	price_.resize(mat_->getNumRows());
	activity_.resize(mat_->getNumRows());

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
		rlbd_.push_back(rlbd0[j]);
		rubd_.push_back(rubd0[j]);
		sense_.push_back(sense0[j]);
		rhs_.push_back(rhs0[j]);
		range_.push_back(range0[j]);
	}

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
	clbd_.assign(collb, collb + mat_->getNumCols());
	cubd_.assign(colub, colub + mat_->getNumCols());
	obj_.assign(obj, obj + mat_->getNumCols());
	rlbd_.assign(rowlb, rowlb + mat_->getNumRows());
	rubd_.assign(rowub, rowub + mat_->getNumRows());
	sense_.resize(mat_->getNumRows());
	rhs_.resize(mat_->getNumRows());
	range_.resize(mat_->getNumRows());
	for (int i = 0; i < mat_->getNumRows(); ++i) {
		DSPdebugMessage("Load row %d: rlbd [%e] rubd [%e]\n", i, rowlb[i], rowub[i]);
		convertBoundToSense(rowlb[i], rowub[i], sense_[i], rhs_[i], range_[i]);
	}
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
	clbd_.assign(collb, collb + mat_->getNumCols());
	cubd_.assign(colub, colub + mat_->getNumCols());
	obj_.assign(obj, obj + mat_->getNumCols());
	sense_.assign(rowsen, rowsen + mat_->getNumRows());
	rhs_.assign(rowrhs, rowrhs + mat_->getNumRows());
	range_.assign(rowrng, rowrng + mat_->getNumRows());
	rlbd_.resize(mat_->getNumRows());
	rubd_.resize(mat_->getNumRows());
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
	status_(UNKNOWN),
	objval_(0.0),
	nIters_(0),
	qpgen_(NULL),
	prob_(NULL),
	vars_(NULL),
	resid_(NULL),
	solver_(NULL),
	nnzQ_(0),
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
	loadProblem(*(rhs.mat_), &rhs.clbd_[0], &rhs.cubd_[0], &rhs.obj_[0],
			rhs.nnzQ_, &rhs.irowQ_[0], &rhs.jcolQ_[0], &rhs.dQ_[0], &rhs.rlbd_[0], &rhs.rubd_[0]);
	/** copy solution data */
	x_ = rhs.x_;
	reduced_ = rhs.reduced_;
	price_ = rhs.price_;
	activity_ = rhs.activity_;
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
		loadProblem(*(rhs.mat_), &rhs.clbd_[0], &rhs.cubd_[0], &rhs.obj_[0],
				rhs.nnzQ_, &rhs.irowQ_[0], &rhs.jcolQ_[0], &rhs.dQ_[0], &rhs.rlbd_[0], &rhs.rubd_[0]);
		/** copy solution data */
		x_ = rhs.x_;
		reduced_ = rhs.reduced_;
		price_ = rhs.price_;
		activity_ = rhs.activity_;
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
