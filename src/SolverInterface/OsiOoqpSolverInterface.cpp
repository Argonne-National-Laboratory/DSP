/*
 * OsiOoqpSolverInterface.cpp
 *
 *  Created on: Sep 9, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include <memory>

#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"

#include "Status.h"
#include "QpGenSparseMa27.h"
#include "GondzioSolver.h"
#include "MehrotraSolver.h"
#include "SimpleVector.h"

#include "Utility/DspMessage.h"
#include "SolverInterface/OsiOoqpSolverInterface.hpp"

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
	int my = 0;   /**< # of equality constraints */
	int mz = 0;   /**< # of inequality constraints */
	int nnzA = 0; /**< # of nonzero elements of the equality constraints */
	int nnzC = 0; /**< # of nonzero elements of the inequality constraints */

	/** row-wise matrix */
	const CoinPackedMatrix* mat = getMatrixByRow();
	DSPdebug(mat->verifyMtx(4));

	/**
	 * The duplicate rows of the matrix will be removed,
	 * as OOQP may result in numerical issues with the duplicates.
	 * The resulting matrix will be stored in reduced_mat with mapping rows.
	 */

	/** create a matrix after reduction */
	std::shared_ptr<CoinPackedMatrix> reduced_mat(new CoinPackedMatrix(false, 0, 0));
	reduced_mat->setDimensions(0, nx);

	/** map the reduced matrix row index to the original matrix row index */
	std::map<int,int> redmat_to_orglow, redmat_to_orgupp;
	/** row bounds of the reduced matrix */
	std::vector<double> reduced_rlbd, reduced_rubd;

	/** check duplicate rows; and create the reduced matrix */
	for (int i = 0; i < mat->getNumRows(); ++i) {
		bool dup = false;
		const CoinShallowPackedVector row = mat->getVector(i);
		if (row.getNumElements() == 0) continue;
#ifdef DSP_DEBUG
		printf("row %d: ", i);
		for (int j = 0; j < row.getNumElements(); ++j)
			printf("[%d] %+e ", row.getIndices()[j], row.getElements()[j]);
		printf("rlbd %+e rubd %+e\n", rlbd_[i], rubd_[i]);
#endif
		for (int ired = 0; ired < reduced_mat->getNumRows(); ++ired) {
			const CoinShallowPackedVector urow = reduced_mat->getVector(ired);
			double scale_factor = 1.0;

			/** found a duplicate? */
			dup = false;
			if (urow.getNumElements() == row.getNumElements() &&
					std::equal(urow.getIndices(),
							urow.getIndices() + urow.getNumElements(),
							row.getIndices())) {
				if (std::equal(urow.getElements(),
							urow.getElements() + urow.getNumElements(),
							row.getElements()))
					dup = true;
				else {
					/** Need scaling */
					scale_factor = urow.getElements()[0] / row.getElements()[0];
					dup = true;
					for (int j = 0; j < urow.getNumElements(); ++j) {
						double coef_diff = urow.getElements()[j] - scale_factor * row.getElements()[j];
						if (fabs(coef_diff) > 1.0e-8) {
							dup = false;
							break;
						}
					}
				}
			}

			if (dup) {
				double scaled_rlbd, scaled_rubd;
				if (scale_factor > 0) {
					scaled_rlbd = scale_factor * rlbd_[i];
					scaled_rubd = scale_factor * rubd_[i];
				} else {
					scaled_rlbd = scale_factor * rubd_[i];
					scaled_rubd = scale_factor * rlbd_[i];
				}

				if (reduced_rlbd[ired] < reduced_rubd[ired]) {
					/** found a tighter lower bound? */
					if (reduced_rlbd[ired] < scaled_rlbd) {
						reduced_rlbd[ired] = scaled_rlbd;
						redmat_to_orglow[ired] = i;
					}
					/** found a tighter upper bound? */
					if (reduced_rubd[ired] > scaled_rubd) {
						reduced_rubd[ired] = scaled_rubd;
						redmat_to_orgupp[ired] = i;
					}
				} else if (reduced_rlbd[ired] >= scaled_rlbd - 1.0e-8 &&
						reduced_rlbd[ired] <= scaled_rubd + 1.0e-8) {
					/** nothing to do */
				} else {
					/** force to add and trigger infeasibility */
					dup = false;
				}
			}

			if (dup) break;
		}
		/** found a unique row? */
		if (!dup) {
			reduced_mat->appendRow(row.getNumElements(), row.getIndices(), row.getElements());
			reduced_rlbd.push_back(rlbd_[i]);
			reduced_rubd.push_back(rubd_[i]);
			redmat_to_orglow[reduced_mat->getNumRows()-1] = i;
			redmat_to_orgupp[reduced_mat->getNumRows()-1] = i;
		}
	}

	/** count nonzeros; and create the row maps */
	ia_to_orglow_.clear();
	ia_to_orgupp_.clear();
	iclow_to_org_.clear();
	icupp_to_org_.clear();
	for (int i = 0; i < reduced_mat->getNumRows(); ++i) {
		const CoinShallowPackedVector row = reduced_mat->getVector(i);
		if (reduced_rlbd[i] == reduced_rubd[i]) {
			ia_to_orglow_.push_back(redmat_to_orglow[i]);
			ia_to_orgupp_.push_back(redmat_to_orgupp[i]);
			my++;
			nnzA += row.getNumElements();
		} else {
			iclow_to_org_.push_back(redmat_to_orglow[i]);
			icupp_to_org_.push_back(redmat_to_orgupp[i]);
			mz++;
			nnzC += row.getNumElements();
		}
	}
#ifdef DSP_DEBUG
	for (int i = 0; i < reduced_mat->getNumRows(); ++i)
		printf("i [%d]: redmat_to_orglow %d, redmat_to_orgupp %d, reduced_rlbd %+e, reduced_rubd %+e, # of elements %d\n",
				i, redmat_to_orglow[i], redmat_to_orgupp[i], reduced_rlbd[i], reduced_rubd[i], reduced_mat->getVector(i).getNumElements());

	for (int i = 0; i < reduced_mat->getNumRows(); ++i) {
		const CoinShallowPackedVector row = reduced_mat->getVector(i);
		printf("row %d: ", i);
		for (int j = 0; j < row.getNumElements(); ++j)
			printf("[%d] %+e ", row.getIndices()[j], row.getElements()[j]);
		printf("rlbd %+e rubd %+e\n", reduced_rlbd[i], reduced_rubd[i]);
	}
#endif

	DSPdebugMessage("nx(# of cols) %d, my(# of Eq.) %d, mz (# of Ineq.) %d, nnzA %d, nnzC %d, unique rows %d\n",
			nx, my, mz, nnzA, nnzC, reduced_mat->getNumRows());

	std::vector<double> xlow(nx, 0.0);
	std::vector<char> ixlow(nx, 0);
	std::vector<double> xupp(nx, 0.0);
	std::vector<char> ixupp(nx, 0);
	std::vector<int> irowA;
	std::vector<int> jcolA;
	std::vector<double> dA;
	std::vector<double> b(my, 0.0);
	std::vector<int> irowC;
	std::vector<int> jcolC;
	std::vector<double> dC;
	std::vector<double> clow(mz, 0.0);
	std::vector<char> iclow(mz, 0);
	std::vector<double> cupp(mz, 0.0);
	std::vector<char> icupp(mz, 0);

	/** reserve memory */
	irowA.reserve(nnzA);
	jcolA.reserve(nnzA);
	dA.reserve(nnzA);
	irowC.reserve(nnzC);
	jcolC.reserve(nnzC);
	dC.reserve(nnzC);

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
	}

	int iposA = 0, iposC = 0;
	for (int i = 0; i < reduced_mat->getNumRows(); ++i) {
		const CoinShallowPackedVector row = reduced_mat->getVector(i);
		if (reduced_rlbd[i] == reduced_rubd[i]) {
			for (int j = 0; j < row.getNumElements(); ++j) {
				irowA.push_back(iposA);
				jcolA.push_back(row.getIndices()[j]);
				dA.push_back(row.getElements()[j]);
#ifdef DSP_DEBUG_MORE
				printf("i %d, irowA %d, jcolA %d, dA %e\n", i, iposA, row.getIndices()[j], row.getElements()[j]);
#endif
			}
			b[iposA] = reduced_rlbd[i];
			iposA++;
#ifdef DSP_DEBUG_MORE
			printf("b %e\n", reduced_rlbd[i]);
#endif
		} else {
			for (int j = 0; j < row.getNumElements(); ++j) {
				irowC.push_back(iposC);
				jcolC.push_back(row.getIndices()[j]);
				dC.push_back(row.getElements()[j]);
			}
			if (reduced_rubd[i] < inf) {
				cupp[iposC] = reduced_rubd[i];
				icupp[iposC] = 1;
			}
			if (reduced_rlbd[i] > -inf) {
				clow[iposC] = reduced_rlbd[i];
				iclow[iposC] = 1;
			}
			iposC++;
		}
	}

	/** create QP generator */
	qpgen = new QpGenSparseMa27(nx, my, mz, nnzQ_, dA.size(), dC.size());
	//qpgen = new QpGenSparseMa57(nx, my, mz, nnzQ_, nnzA, nnzC);

	std::vector<double> obj_sensed(obj_);
	for (int j = 0; j < nx; ++j)
		obj_sensed[j] *= obj_sense_;

	prob = (QpGenData*)dynamic_cast<QpGenSparseSeq*>(qpgen)->copyDataFromSparseTriple(
			&obj_sensed[0],  &irowQ_[0], nnzQ_,     &jcolQ_[0], &dQ_[0],
			&xlow[0],  &ixlow[0],  &xupp[0],  &ixupp[0],
			&irowA[0], nnzA,       &jcolA[0], &dA[0],     &b[0],
			&irowC[0], nnzC,       &jcolC[0], &dC[0],
			&clow[0],  &iclow[0],  &cupp[0],  &icupp[0]);
	DSPdebug(prob->print());
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
	//solver_ = new GondzioSolver(qpgen_, prob_);
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

	/** FIXME: Why doesn't this work? */
	if (messageHandler()->logLevel() >= 4)
		solver_->monitorSelf();

	/** solve */
	status_ = solver_->solve(prob_, vars_, resid_);
	DSPdebugMessage("logLevel %d, OOQP status %d\n", messageHandler()->logLevel(), status_);

	if (status_ == SUCCESSFUL_TERMINATION || status_ == MAX_ITS_EXCEEDED) {

		nIters_ = solver_->iter;

		/** allocate memory for Osi solutions */
		freeCachedResults();
		x_.resize(prob_->nx);
		reduced_.resize(prob_->nx, 0.0);
		activity_.resize(mat_->getNumRows(), 0.0);
		price_.resize(mat_->getNumRows(), 0.0);

		/** retrieve OOQP solutions */
		objval_ = prob_->objectiveValue(vars_);
		vars_->x->copyIntoArray(&x_[0]);

		double* iclow  = dynamic_cast<SimpleVector*>(prob_->iclow.ptr())->elements();
		double* icupp  = dynamic_cast<SimpleVector*>(prob_->icupp.ptr())->elements();
		double* ixlow  = dynamic_cast<SimpleVector*>(prob_->ixlow.ptr())->elements();
		double* ixupp  = dynamic_cast<SimpleVector*>(prob_->ixupp.ptr())->elements();
		double* y      = dynamic_cast<SimpleVector*>(vars_->y.ptr())->elements();      /**< dual variables corresponding to equality constraints */
		double* lambda = dynamic_cast<SimpleVector*>(vars_->lambda.ptr())->elements(); /**< dual variables corresponding to inequality (>=) constraints */
		double* pi     = dynamic_cast<SimpleVector*>(vars_->pi.ptr())->elements();     /**< dual variables corresponding to inequality (<=) constraints */
		double* gamma  = dynamic_cast<SimpleVector*>(vars_->gamma.ptr())->elements();  /**< dual variables corresponding to column lower bounds */
		double* phi    = dynamic_cast<SimpleVector*>(vars_->phi.ptr())->elements();    /**< dual variables corresponding to column upper bounds */

		/** calculate activity */
		mat_->times(&x_[0], &activity_[0]);

		/** calculate price */
		for (int i = 0; i < prob_->my; ++i) {
			//printf("y[%d] %+e: ia_to_orglow_ %d, ia_to_orgupp_ %d\n", i, y[i], ia_to_orglow_[i], ia_to_orgupp_[i]);
			if (y[i] > 0)
				price_[ia_to_orglow_[i]] += y[i];
			else
				price_[ia_to_orgupp_[i]] += y[i];
		}
		for (int i = 0; i < prob_->mz; ++i) {
			//printf("mz[%d]: lambda %+e, pi %+e, iclow_to_org_ %d, icupp_to_org_ %d\n", i, lambda[i], pi[i], iclow_to_org_[i], icupp_to_org_[i]);
			if (iclow[i] > 0)
				price_[iclow_to_org_[i]] += lambda[i];
			if (icupp[i] > 0)
				price_[icupp_to_org_[i]] -= pi[i];
		}

		for (int j = 0; j < prob_->nx; ++j) {
			if (ixlow[j] > 0)
				reduced_[j] += gamma[j];
			if (ixupp[j] > 0)
				reduced_[j] -= phi[j];
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
	return objval_ * obj_sense_;
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

	if (obj_.size() != mat_->getNumCols())
		throw CoinError("obj_ size is not compatible to the number of columns in the matrix", "deleteCols", "OsiOoqpSolverInterface");

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
#ifdef DSP_DEBUG_MORE
		printf("Load row %d: rlbd [%e] rubd [%e]\n", i, rowlb[i], rowub[i]);
#endif
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
