/*
 * PipsInterface.h
 *  Created on: Dec 12, 2017
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_PIPSINTERFACE_H_
#define SRC_SOLVER_DANTZIGWOLFE_PIPSINTERFACE_H_

#include "stochasticInput.hpp"
#include "CoinFinite.hpp"
#include "CoinPackedVector.hpp"
//#include "Solver/DantzigWolfe/DwBundleDualPips.h"
//#include "Solver/DantzigWolfe/DwWorkerPips.h"

class DwBundleDualPips;
class DwWorkerPips;

class PipsInput: public stochasticInput {
	friend class PipsInterface;
private:
	int nscen_;
	int nvars1_;
	int ncons1_;
	int nvars2_;
	std::vector<int> ncons2_;
	std::vector<double> clbd1_;
	std::vector<double> cubd1_;
	std::vector<double> obj1_;
	std::vector<std::string> cname1_;
	std::vector<double> rlbd1_;
	std::vector<double> rubd1_;
	std::vector<std::string> rname1_;
	std::vector<double> clbd2_;
	std::vector<double> cubd2_;
	std::vector<double> obj2_;
	std::vector<std::vector<std::string> > cname2_;
	std::vector<std::vector<double> > rlbd2_;
	std::vector<std::vector<double> > rubd2_;
	std::vector<std::vector<std::string> > rname2_;
	CoinPackedMatrix matA_;
	std::vector<CoinPackedMatrix> matW_;
	std::vector<CoinPackedMatrix> matT_;
	CoinPackedMatrix hess1_;

public:
	/** default constructor */
	PipsInput(int nscen, int nvars1) : nscen_(nscen), nvars1_(nvars1), ncons1_(0) {
		nvars2_ = 1 + nvars1_;
		ncons2_.resize(nscen_, 0);
		clbd1_.resize(nvars1_, -COIN_DBL_MAX);
		cubd1_.resize(nvars1_, +COIN_DBL_MAX);
		obj1_.resize(nvars1_, 0.0);
		for (int j = 0; j < nvars1_; ++j)
			cname1_.push_back("c" + std::to_string(j));
		for (int j = 0; j < ncons1_; ++j)
			rname1_.push_back("r" + std::to_string(j));
		clbd2_.resize(nvars2_, -COIN_DBL_MAX);
		cubd2_.resize(nvars2_, +COIN_DBL_MAX);
		obj2_.resize(nvars2_, 0.0);
		obj2_[0] = -1.0;
		rlbd2_.resize(nscen_);
		rubd2_.resize(nscen_);
		cname2_.resize(nscen_); rname2_.resize(nscen_);
		matA_.setDimensions(ncons1_, nvars1_);
		matW_.resize(nscen_, CoinPackedMatrix(false,0,0));
		matT_.resize(nscen_, CoinPackedMatrix(false,0,0));
		for (int j = 0; j < nscen_; ++j) {
			matW_[j].setDimensions(ncons2_[j], nvars2_);
			matT_[j].setDimensions(ncons2_[j], nvars1_);
			for (int i = 0; i < nvars1_; ++i) {
				int cind = i; double cval = 1.0;
				matT_[j].appendRow(1, &cind, &cval);
				cind = i+1; cval = -1.0;
				matW_[j].appendRow(1, &cind, &cval);
				rlbd2_[j].push_back(0.0);
				rubd2_[j].push_back(0.0);
				ncons2_[j]++;
			}
			for (int k = 0; k < nvars2_; ++k)
				cname2_[j].push_back("s" + std::to_string(j) + "c" + std::to_string(k));
			for (int k = 0; k < ncons2_[j]; ++k)
				rname2_[j].push_back("s" + std::to_string(j) + "r" + std::to_string(k));
			//matW_[j].verifyMtx(4);
			//matT_[j].verifyMtx(4);
		}
		hess1_.setDimensions(nvars1_, nvars1_);
	}

	/** copy constructor */
	PipsInput(const PipsInput& rhs):
		nscen_(rhs.nscen_), nvars1_(rhs.nvars1_), ncons1_(rhs.ncons1_), nvars2_(rhs.nvars2_), ncons2_(rhs.ncons2_), 
		obj1_(rhs.obj1_), cname1_(rhs.cname1_), rname1_(rhs.rname1_), 
		clbd2_(rhs.clbd2_), cubd2_(rhs.cubd2_), obj2_(rhs.obj2_), cname2_(rhs.cname2_),
		rlbd2_(rhs.rlbd2_), rubd2_(rhs.rubd2_), rname2_(rhs.rname2_),
		matA_(rhs.matA_), matW_(rhs.matW_), matT_(rhs.matT_),
		hess1_(rhs.hess1_) {} 
		
	//virtual ~PipsInput() {}
	virtual int nScenarios() { return nscen_; }
	virtual int nFirstStageVars() { return nvars1_; }
	virtual int nFirstStageCons() { return ncons1_; }
	virtual int nSecondStageVars(int scen) { return nvars2_; }
	virtual int nSecondStageCons(int scen) { return ncons2_[scen]; }

	virtual std::vector<double> getFirstStageColLB() { return clbd1_; }
	virtual std::vector<double> getFirstStageColUB() { return cubd1_; }
	virtual std::vector<double> getFirstStageObj() { return obj1_; }
	virtual std::vector<std::string> getFirstStageColNames() { return cname1_; }
	virtual std::vector<double> getFirstStageRowLB() { return rlbd1_; }
	virtual std::vector<double> getFirstStageRowUB() { return rubd1_; }
        //virtual std::vector<double> getLinkRowLB(){ return std::vector<double>(); } /** TODO */
        //virtual std::vector<double> getLinkRowUB(){ return std::vector<double>(); } /** TODO */
	virtual std::vector<std::string> getFirstStageRowNames() { return rname1_; }
	virtual bool isFirstStageColInteger(int col) { return false; } /** why do we need this? */

	virtual std::vector<double> getSecondStageColLB(int scen) { return clbd2_; }
	virtual std::vector<double> getSecondStageColUB(int scen) { return cubd2_; }
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen) { return obj2_; }
	virtual std::vector<std::string> getSecondStageColNames(int scen) { return cname2_[scen]; }
	virtual std::vector<double> getSecondStageRowUB(int scen) { return rubd2_[scen]; }
	virtual std::vector<double> getSecondStageRowLB(int scen) { return rlbd2_[scen]; }
	virtual std::vector<std::string> getSecondStageRowNames(int scen) { return rname2_[scen]; }
	virtual double scenarioProbability(int scen) { return (1.0/nscen_); }
	virtual bool isSecondStageColInteger(int scen, int col) { return false; } /** why do we need this? */

	// returns the column-oriented first-stage constraint matrix (A matrix) 
	virtual CoinPackedMatrix getFirstStageConstraints() { return matA_; }
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) { return matW_[scen]; }
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen) { return matT_[scen]; }

	// some problem characteristics that could be helpful to know
	
	// all scenarios have the same number of variables and constraints
	virtual bool scenarioDimensionsEqual() { return false; } /** TODO: Can this be false? */
	// constraint (and hessian) matrices are identical for each scenario,
	// column and row bounds and objective are allowed to vary
	virtual bool onlyBoundsVary() { return false; }
	// all scenarios equally likely
	virtual bool allProbabilitiesEqual() { return true; }
	// all second-stage variables continuous
	virtual bool continuousRecourse() { return true; }

	/* Quadratic terms:
	We allow the input to specify quadratic objective in the form of:
	(1/2)x^TQx + c^T + \sum_{i=1}^N p_i ((1/2)y_i^TQ_iy_i + x_i^T{\hat Q_i^T}y_i + c_i^Ty)
	Q is the first-stage hessian
	Q_i is the second-stage hessian
	\hat Q_i is the second-stage cross hessian
	Default implementations are provided so that these do not need to be implemented if not used
	*/

	// column-oriented *lower triangle only*
	// Q
	virtual CoinPackedMatrix getFirstStageHessian() { return hess1_; }

	/* TODO
        virtual int nLinkCons() { return 0; } 
        virtual int nLinkECons() { return 0; }
        virtual int nLinkICons() { return 0; }
        virtual CoinPackedMatrix getLinkMatrix(int nodeid){ return CoinPackedMatrix(); }
	std::string datarootname;
	int useInputDate;
	*/
};

class PipsInterface {
	friend class DwBundleDualPips;
	friend class DwWorkerPips;
private:
	PipsInput* input_;
	int nrows_; /**< number of rows */

	/** objective value and solutions */
	double objval_;
	std::vector<double> solution_;

public:
	std::vector<int> scenarios_;
	std::vector<double> thetas_;

public:
	/** constructor */
	PipsInterface(int nscen, int ncols);

	/** copy constructor */
	PipsInterface(const PipsInterface& rhs): objval_(rhs.objval_), solution_(rhs.solution_), nrows_(rhs.nrows_) {
		input_ = new PipsInput(*(rhs.input_));
	}

	/** destructor */
	virtual ~PipsInterface() {
		delete input_;
	}

	virtual int getNumRows() { return nrows_; }
	virtual int getNumCols() { return (input_->nvars1_ + input_->nscen_ * input_->nvars2_); }
	virtual int getNumScenarios() { return input_->nscen_; }

	virtual void setObjCoef(int j, double val) {
		input_->obj1_[j] = val;
	}

	/** retrieve Bundle data to distribute */
	virtual void retrieveBundles(std::vector<int>& rstart, std::vector<CoinPackedVector*>& vecs, std::vector<double>& rlbd, std::vector<double>& rubd, std::vector<int>& scen);

	/** complete PipsInterface data and solve */
	virtual int solve(double weight);

	virtual void addRow(const CoinPackedVector& v, double lb, double ub);

	/** remove Bundle cuts for matrices T and W */
	virtual void clearMatricesTW();

	/** print all data */
	virtual void print();
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_PIPSINTERFACE_H_ */
