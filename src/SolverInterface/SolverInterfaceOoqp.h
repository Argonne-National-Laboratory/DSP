/*
 * SolverInterfaceOoqp.h
 *
 *  Created on: Jan 7, 2015
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_SOLVERINTERFACEOOQP_H_
#define SRC_SOLVER_SOLVERINTERFACEOOQP_H_

//#define DSP_HAS_MA57

/** DSP */
#include <Utility/DspMessage.h>
#include "SolverInterface/SolverInterface.h"
#include "QpGen.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "Solver.h"

#define GUTS_OF_LOAD_PROBLEM gutsOfLoadSparseProblem
#define SET_HESSIAN          setSparseHessian

/** columns dynamically added */
class dynColumns
{
public:
	dynColumns(int nrows, int ncols)
	{
		cols_ = new CoinPackedMatrix();
		cols_->setDimensions(nrows, 0);
	}

	dynColumns(dynColumns * dynCols)
	{
		assert(dynCols);
		if (dynCols->cols_->getNumCols() > 0)
		{
			cols_ = new CoinPackedMatrix(*(dynCols->cols_));
			objs_ = dynCols->objs_;
			lbs_ = dynCols->lbs_;
			ubs_ = dynCols->ubs_;
		}
		else
		{
			cols_ = new CoinPackedMatrix();
			cols_->setDimensions(dynCols->cols_->getNumRows(), 0);
		}
	}

	~dynColumns()
	{
		FREE_PTR(cols_);
	}

	void addCol(int size, const int * indices, const double * vals, double lb, double ub, double obj)
	{
		objs_.push_back(obj);
		lbs_.push_back(lb);
		ubs_.push_back(ub);
		cols_->appendCol(size, indices, vals);
	}

	void clear()
	{
		objs_.clear();
		lbs_.clear();
		ubs_.clear();
		cols_->deleteMajorVectors(cols_->getMajorDim(), cols_->getMajorIndices());
	}

	int size() {return cols_->getNumCols();}

	vector<double> objs_;     /**< objective coefficient */
	CoinPackedMatrix * cols_; /**< column */
	vector<double> lbs_;      /**< column lower bound */
	vector<double> ubs_;      /**< column upper bound */
};

class SolverInterfaceOoqp: public SolverInterface
{
public:
	SolverInterfaceOoqp(DspParams * par) :
		SolverInterface(par),
		dynCols_(NULL),
		sense_(1),
		status_(DSP_STAT_UNKNOWN),
		print_level_(0),
		qp_(NULL),
		prob_(NULL),
		vars_(NULL),
		resid_(NULL),
		solver_(NULL),
		released_(true),
		nx_(0), my_(0), mz_(0),
		nrows_(0), ncols_(0),
		xlow_(NULL), ixlow_(NULL), xupp_(NULL), ixupp_(NULL),
		c_(NULL),
		mat_(NULL),
		clbd_(NULL), cubd_(NULL),
		obj_(NULL),
		rlbd_(NULL), rubd_(NULL),
		nnzA_(0), irowA_(NULL), jcolA_(NULL), dA_(NULL), b_(NULL),
		nnzC_(0), irowC_(NULL), jcolC_(NULL), dC_(NULL),
		clow_(NULL), iclow_(NULL), cupp_(NULL), icupp_(NULL),
		nnzQ_(0), irowQ_(NULL), jcolQ_(NULL), dQ_(NULL),
		objval_(COIN_DBL_MAX),
		x_(NULL),
		y_(NULL),
		lambda_(NULL),
		pi_(NULL),
		gamma_(NULL),
		phi_(NULL),
		nLambdas_(0),
		nPis_(0),
		nIters_(0) {}

	/** copy constructor */
	SolverInterfaceOoqp(SolverInterfaceOoqp * si);

	/** clone */
	virtual SolverInterface * clone();

	virtual ~SolverInterfaceOoqp();

	/** load problem */
	virtual void loadProblem(
			OsiSolverInterface * si,
			const char * probname = "null");

	/** load problem */
	virtual void loadProblem(
			CoinPackedMatrix * mat,
			const double * collb,
			const double * colub,
			const double * obj,
			const char * ctype,
			const double * rowlb,
			const double * rowub,
			const char * probname = "null");

	/** core part for load problem */
	virtual void gutsOfLoadSparseProblem();

	/** core part for load problem */
	virtual void gutsOfLoadDenseProblem();

	/** add row */
	virtual void addRow(int size, const int * indices, const double * vals, double lb, double ub);

	/** add column */
	virtual void addCol(int size, const int * indices, const double * vals, double lb, double ub, double obj);

	/** clear columns */
	virtual void clearCols();

	/** add cuts */
	virtual int addCuts(OsiCuts cuts, double effectiveness = 0.0);

	/** solve */
	virtual void solve();

	/** get objective sense */
	virtual int getObjSense() {return sense_;}

	/** solution status */
	virtual DSP_RTN_CODE getStatus() {return status_;}

	/** get number of rows */
	virtual int getNumRows() {return my_ + mz_;}

	/** get number of equality constraints */
	virtual int getNumEqConss() {return my_;}

	/** get number of inequality constraints */
	virtual int getNumIneqConss() {return mz_;}

	/** get number of columns */
	virtual int getNumCols() {return ncols_;}

	/** get number of integer variables */
	virtual int getNumIntegers() {return 0;}

	/** get constraint matrix */
	virtual CoinPackedMatrix * getMatrix() {return mat_;}

	/** get column lower bounds */
	virtual const double * getColLower() {return clbd_;}

	/** get column upper bounds */
	virtual const double * getColUpper() {return cubd_;}

	/** get objective */
	virtual const double * getObjCoef() {return obj_;}

	/** get row lower bounds */
	virtual const double * getRowLower() {return rlbd_;}

	/** get row upper bounds */
	virtual const double * getRowUpper() {return rubd_;}

	/** get global primal bound (upper bound in minimization) */
	virtual double getPrimalBound() {return objval_ * sense_;}

	/** get dual bound (lower bound in minimization) */
	virtual double getDualBound() {return objval_ * sense_;}

	/** get solution values */
	virtual const double * getSolution() {return x_;}

	/** get solution values */
	virtual const double * y() {return y_;}

	/** get number of lambdas */
	virtual int getNumLambdas() {return nLambdas_;}

	/** get number of pis */
	virtual int getNumPis() {return nPis_;}

	/** get lambda */
	virtual const double * lambda() {return lambda_;}

	/** get pi */
	virtual const double * pi() {return pi_;}

	/** get gamma */
	virtual const double * gamma() {return gamma_;}

	/** get phi */
	virtual const double * phi() {return phi_;}

	/** get number of iterations */
	virtual int getIterationCount() {return nIters_;}

	/** get number of nodes */
	virtual int getNumNodes() {return 0;}

	/** get print out level */
	virtual int getPrintLevel() {return print_level_;}

	/** OOQP objects release? */
	virtual bool released() {return released_;}

	/** set objective coefficients */
	virtual void setObjSense(int sense) {sense_ = sense;}

	/** set objective coefficients */
	virtual void setObjCoef(double * obj);

	/** set lower triangular Hessian matrix */
	virtual void setHessian(int nnzQ, int * irowQ, int * jcolQ, double * dQ);

	/** set column bounds */
	virtual void setColLower(int index, double lb);

	/** set column bounds */
	virtual void setColUpper(int index, double ub);

	/** set column bounds */
	virtual void setColBounds(int index, double lb, double ub);

	/** set row bounds */
	virtual void setRowLower(int index, double lb);

	/** set row bounds */
	virtual void setRowUpper(int index, double ub);

	/** set row bounds */
	virtual void setRowBounds(int index, double lb, double ub);

	/** set node limit */
	virtual void setNodeLimit(int limit) {}

	/** set iteration limit */
	virtual void setIterLimit(int limit) {}

	/** set wall time limit */
	virtual void setTimeLimit(double sec) {DSPdebugMessage("Not available for OOQP.\n");}

	/** set gap tolerance */
	virtual void setGapTol(double tol) {}

	/** set print out level */
	virtual void setPrintLevel(int level) {print_level_ = level;}

	/** write MPS file */
	virtual void writeMps(const char * filename) {DSPdebugMessage("Not available for OOQP.\n");}

protected:

	/** set lower triangular Hessian matrix */
	virtual void setSparseHessian(int nnzQ, int * irowQ, int * jcolQ, double * dQ);

	/** set lower triangular Hessian matrix */
	virtual void setDenseHessian(int nnzQ, int * irowQ, int * jcolQ, double * dQ);

	/** initialize solver interface */
	virtual DSP_RTN_CODE initialize() {return DSP_RTN_OK;}

	/** finalize solver interface */
	virtual DSP_RTN_CODE finalize();

	/** release OOQP objects */
	virtual void releaseOOQP();

	/** release solver */
	void release();

public:

	OsiCuts cuts_; /**< for additional rows/columns */

	dynColumns * dynCols_; /**< dynamically added columns */

protected:

	int sense_; /**< objective sense */
	DSP_RTN_CODE status_; /**< solution status */
	int print_level_;

	/** OOQP objects */
	QpGen * qp_;
	QpGenData * prob_;
	QpGenVars * vars_;
	QpGenResiduals * resid_;
	Solver * solver_;
	bool released_;

	/** original data */
	int nx_;
	int my_;
	int mz_;
	int nrows_;
	int ncols_;
	double * xlow_;
	char * ixlow_;
	double * xupp_;
	char * ixupp_;
	double * c_;
	CoinPackedMatrix * mat_;
	double * clbd_;
	double * cubd_;
	double * obj_;
	double * rlbd_;
	double * rubd_;
	int nnzA_;
	int * irowA_;
	int * jcolA_;
	double * dA_;
	double * b_;
	int nnzC_;
	int * irowC_;
	int * jcolC_;
	double * dC_;
	double * clow_;
	char * iclow_;
	double * cupp_;
	char * icupp_;
	int nnzQ_;
	int * irowQ_;
	int * jcolQ_;
	double * dQ_;

	double objval_;   /**< objective value */
	double * x_;      /**< primal solution */
	double * y_;      /**< dual variables corresponding to equality constraints */
	double * lambda_; /**< dual variables corresponding to inequality (>=) constraints */
	double * pi_;     /**< dual variables corresponding to inequality (<=) constraints */
	double * gamma_;  /**< dual variables corresponding to column lower bounds */
	double * phi_;    /**< dual variables corresponding to column upper bounds */
	int nLambdas_;    /**< number of lambdas */
	int nPis_;        /**< number of pis */
	int nIters_;      /**< number of iterations */
};

#endif /* SRC_SOLVER_SOLVERINTERFACEOOQP_H_ */
