/*
 * PipsInterface.h
 *  Created on: Dec 12, 2017
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_PIPSINTERFACE_H_
#define SRC_SOLVER_DANTZIGWOLFE_PIPSINTERFACE_H_

#include <numeric>
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
	std::vector<int> nvars2_;
	int ncons2_;
	std::vector<double> clbd1_;
	std::vector<double> cubd1_;
	std::vector<double> obj1_;
	std::vector<std::string> cname1_;
	std::vector<double> rlbd1_;
	std::vector<double> rubd1_;
	std::vector<std::string> rname1_;
	std::vector<std::vector<double>> clbd2_;
	std::vector<std::vector<double>> cubd2_;
	std::vector<std::vector<double>> obj2_;
	std::vector<std::vector<double>> cx2_;
	std::vector<std::vector<CoinPackedVector>> x2_;
	std::vector<std::vector<std::string> > cname2_;
	std::vector<double> rlbd2_;
	std::vector<double> rubd2_;
	std::vector<std::vector<std::string> > rname2_;
	CoinPackedMatrix matA_;
	CoinPackedMatrix matT_;
	std::vector<CoinPackedMatrix> matW_;
	CoinPackedMatrix* hessxx_;
	std::vector<CoinPackedMatrix*> hessxy_; /**< Hessian associated with both first- and second-stage variables */
	std::vector<CoinPackedMatrix*> hessyy_; /**< Hessian associated with both first- and second-stage variables */

public:
	/** default constructor */
	PipsInput(int nscen) : nscen_(nscen) {

		if (nvars1_ > 0) {
			clbd1_.resize(nvars1_, -COIN_DBL_MAX);
			cubd1_.resize(nvars1_, +COIN_DBL_MAX);
			obj1_.resize(nvars1_, 0.0);
			cname1_.reserve(nvars1_);
			for (int j = 0; j < nvars1_; ++j)
				cname1_.push_back("c" + std::to_string(j));
		}
		ncons2_ = 1;
		rlbd2_.resize(ncons2_, 1.0);
		rubd2_.resize(ncons2_, 1.0);
		nvars2_.resize(nscen_, 0);
		clbd2_.resize(nscen_); cubd2_.resize(nscen_); 
		obj2_.resize(nscen_); cx2_.resize(nscen_); x2_.resize(nscen_);
		cname2_.resize(nscen_); rname2_.resize(nscen_);
		matA_.setDimensions(ncons1_, nvars1_);
		matT_.setDimensions(ncons2_, nvars1_);
		matW_.resize(nscen_, CoinPackedMatrix(true,0,0));
		hessxx_.setDimensions(nvars1_, 0);
		for (int i = 0; i < nvars1_; ++i) {
			int hessi = i;
			double hesse = 1.0 * nscen_;
			hessxx_.appendCol(1, &hessi, &hesse);
		}
		// hessxy_.resize(nscen_, CoinPackedMatrix(false,0,0));
		hessxy_.resize(nscen_, CoinPackedMatrix(true,0,0));
		hessyy_.resize(nscen_, CoinPackedMatrix(true,0,0));
		for (int j = 0; j < nscen_; ++j) {
			clbd2_[j].clear(); cubd2_[j].clear(); obj2_[j].clear(); cname2_[j].clear();
			rname2_[j].reserve(ncons2_);
			for (int k = 0; k < ncons2_; ++k)
				rname2_[j].push_back("s" + std::to_string(j) + "r" + std::to_string(k));
			matW_[j].setDimensions(ncons2_, 0);
			// hessxy_[j].setDimensions(0, nvars1_);
			hessxy_[j].setDimensions(0, nvars1_);
			hessyy_[j].setDimensions(0, 0);
			//matW_[j].verifyMtx(4);
			//hessxy_[j].verifyMtx(4);
		}
	}

	/** copy constructor */
	PipsInput(const PipsInput& rhs):
		nscen_(rhs.nscen_), nvars1_(rhs.nvars1_), ncons1_(rhs.ncons1_), nvars2_(rhs.nvars2_), ncons2_(rhs.ncons2_), 
		obj1_(rhs.obj1_), cname1_(rhs.cname1_), rname1_(rhs.rname1_), 
		clbd2_(rhs.clbd2_), cubd2_(rhs.cubd2_), obj2_(rhs.obj2_), cname2_(rhs.cname2_),
		rlbd2_(rhs.rlbd2_), rubd2_(rhs.rubd2_), rname2_(rhs.rname2_),
		matA_(rhs.matA_), matW_(rhs.matW_),
		hessxx_(rhs.hessxx_), hessxy_(rhs.hessxy_), hessyy_(rhs.hessyy_) {} 
		
	//virtual ~PipsInput() {}
	virtual int nScenarios() { return nscen_; }
	virtual int nFirstStageVars() { return nvars1_; }
	virtual int nFirstStageCons() { return ncons1_; }
	virtual int nSecondStageVars(int scen) { return nvars2_[scen]; }
	virtual int nSecondStageCons(int scen) { return ncons2_; }

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

	virtual std::vector<double> getSecondStageColLB(int scen) { return clbd2_[scen]; }
	virtual std::vector<double> getSecondStageColUB(int scen) { return cubd2_[scen]; }
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen) { return obj2_[scen]; }
	virtual std::vector<std::string> getSecondStageColNames(int scen) { return cname2_[scen]; }
	virtual std::vector<double> getSecondStageRowUB(int scen) { return rubd2_; }
	virtual std::vector<double> getSecondStageRowLB(int scen) { return rlbd2_; }
	virtual std::vector<std::string> getSecondStageRowNames(int scen) { return rname2_[scen]; }
	virtual double scenarioProbability(int scen) { return 1.0/nscen_; }
	virtual bool isSecondStageColInteger(int scen, int col) { return false; } /** why do we need this? */

	// returns the column-oriented first-stage constraint matrix (A matrix) 
	virtual CoinPackedMatrix getFirstStageConstraints() { return matA_; }
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) { return matW_[scen]; }
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen) { return matT_; }

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
	virtual CoinPackedMatrix getFirstStageHessian() { return *hessxx_; }
    virtual CoinPackedMatrix getSecondStageCrossHessian(int scen){ return *(hessxy_[scen]); }
    virtual CoinPackedMatrix getSecondStageHessian(int scen){ return *(hessyy_[scen]); }

	/* TODO
        virtual int nLinkCons() { return 0; } 
        virtual int nLinkECons() { return 0; }
        virtual int nLinkICons() { return 0; }
        virtual CoinPackedMatrix getLinkMatrix(int nodeid){ return CoinPackedMatrix(); }
	std::string datarootname;
	int useInputDate;
	*/

    virtual void loadProblem(int ncons1, int nvars1, int ncons2, std::vector<int> nvars2,
    	CoinPackedMatrix*& hessxx, CoinPackedMatrix**& hessxy, CoinPackedMatrix**& hessyy,
    	std::vector<std::vector<double>> obj2,
    	CoinPackedMatrix*& W);
	virtual void setNumFirstStageCons(int n) {ncons1_=n;}
	virtual void setNumFirstStageVars(int n) {nvars1_=n;}
	virtual void setNumSecondStageCons(int n) {ncons2_=n;}
	virtual void setNumSecondStageVars(std::vector<int> n) {nvars2_=n;}
	virtual void setSecondStageObj(std::vector<std::vector<double>> obj2) {obj2_=obj2;}
	virtual void setFirstStageHessian(CoinPackedMatrix*& H) {
		FREE_PTR(hessxx_);
		hessxx_ = H;
		H = NULL;
	}
	virtual void setSecondStageCrossHessian(CoinPackedMatrix**& H) {
		for (int s = 0; s < nscen_; ++s) {
			FREE_PTR(hessxy_[s]);
			hessxy_[s] = H[s];
			H[s] = NULL;
		}
	}
	virtual void setSecondStageHessian(CoinPackedMatrix**& H) {
		for (int s = 0; s < nscen_; ++s) {
			FREE_PTR(hessyy_[s]);
			hessyy_[s] = H[s];
			H[s] = NULL;
		}
	}
	virtual void setSecondStageConstraints(CoinPackedMatrix**& W) {
		for (int s = 0; s < nscen_; ++s) {
			FREE_PTR(matW_[s]);
			matW_[s] = W[s];
			W[s] = NULL;
		}
	}
	virtual void setSecondStageRowUB(std::vector<double> rubd2) {rubd2_=rubd2;}
	virtual void setSecondStageRowLB(std::vector<double> rlbd2) {rlbd2_=rlbd2;}
};

class PipsInterface {
	friend class DwBundleDualPips;
	friend class DwWorkerPips;
private:
	PipsInput* input_;
	int nrows_; /**< number of rows */
	double penalty_; /**< penalty applied */

	std::vector<int> scenarios_;

	/** objective value and solutions */
	double objval_;
	CoinPackedVector* soln_w_;
	std::vector<CoinPackedVector*> soln_z_;
	std::vector<CoinPackedVector*> soln_theta_;
	std::vector<CoinPackedVector*> soln_lambda_;

	std::vector<double> dualsol_;
	std::vector<double> bestdualsol_;

public:
	/** constructor */
	PipsInterface(int nscen, int ncols);

	/** copy constructor */
	PipsInterface(const PipsInterface& rhs);

	/** destructor */
	virtual ~PipsInterface();

	virtual void clearSolutions();

	virtual int getNumRows() { return nrows_; }
	virtual int getNumCols() { return (input_->nvars1_ + std::accumulate(input_->nvars2_.begin(), input_->nvars2_.end(), 0)); }
	virtual int getNumFirstStageCols() {return input_->nvars1_;}
	virtual int getNumSecondStageCols(int s) {return input_->nvars2_[s];}
	virtual int getNumScenarios() { return input_->nscen_; }

	virtual void updateCenter(double penalty, std::vector<double>& bestdualsol);

	/** collect column data to distribute */
	virtual void collectColumnData(
		double& penalty, /**< current penalty value */
		std::vector<int>& cstart, /**< index for the first new column to read */
		std::vector<int>& scen, /**< scenario index */
		std::vector<CoinPackedVector*>& x, /**< first-stage solution of subproblems */
		std::vector<double>& cx);

	/** complete PipsInterface data and solve */
	virtual int solve(double weight);

	virtual void addCol(int sind, const CoinPackedVector* x, double cx, std::vector<double>& bestdualsol);

	/** remove Bundle cuts for matrices W */
	virtual void clearMatricesW();

	/** print all data */
	virtual void print();
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_PIPSINTERFACE_H_ */
