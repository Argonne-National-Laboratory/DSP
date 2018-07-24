/*
 * DwBundlePipsInput.h
 *  Created on: July 16, 2018
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEPIPSINPUT_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEPIPSINPUT_H_

#include <numeric>
#include "stochasticInput.hpp"
#include "CoinFinite.hpp"
#include "CoinPackedVector.hpp"

class DwBundlePipsInput: public stochasticInput {
public:

	/** NOTE: dynamic part of the input data will be generated as callback */

	/** static part of the input data */
	int nscen_;
	int nvars1_;
	int ncons1_;
	int ncons2_;

	/** proximal parameter */
	double tau_;

	/** proximal center for each scenario */
	std::vector<std::vector<double>> prox_center_;

private:
	/** bundle information */
	std::vector<std::vector<CoinPackedVector>> bundle_x_;
	std::vector<std::vector<double>> bundle_obj_;

public:

	virtual void addBundleInfo(int s, int nelements, int* indices, double* elements, double obj) {
		CoinPackedVector x(nelements, indices, elements);
		bundle_x_[s].push_back(x);
		bundle_obj_[s].push_back(obj);
	}

	virtual void setTau(double tau) { tau_ = tau; }

	virtual void setProxCenter(std::vector<std::vector<double>> c) { prox_center_ = c; }

public:
	/** default constructor */
	DwBundlePipsInput(int nscen, int nvars1);
	virtual int nScenarios() { return nscen_; }
	virtual int nFirstStageVars() { return nvars1_; }
	virtual int nFirstStageCons() { return ncons1_; }
	virtual int nSecondStageVars(int scen);
	virtual int nSecondStageCons(int scen) { return ncons2_; }

	virtual std::vector<double> getFirstStageColLB() { return std::vector<double>(nvars1_, -COIN_DBL_MAX); }
	virtual std::vector<double> getFirstStageColUB() { return std::vector<double>(nvars1_, +COIN_DBL_MAX); }
	virtual std::vector<double> getFirstStageObj() { return std::vector<double>(nvars1_, 0.0); }
	virtual std::vector<std::string> getFirstStageColNames() { return std::vector<std::string>(nvars1_, "w"); }
	virtual std::vector<double> getFirstStageRowLB() { return std::vector<double>(0); }
	virtual std::vector<double> getFirstStageRowUB() { return std::vector<double>(0); }
	virtual std::vector<std::string> getFirstStageRowNames() { return std::vector<std::string>(0); }
	virtual bool isFirstStageColInteger(int col) { return false; }

	virtual std::vector<double> getSecondStageColLB(int scen);
	virtual std::vector<double> getSecondStageColUB(int scen);
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen);
	virtual std::vector<std::string> getSecondStageColNames(int scen);
	virtual std::vector<double> getSecondStageRowUB(int scen) { return std::vector<double>(1, 1.0); }
	virtual std::vector<double> getSecondStageRowLB(int scen) { return std::vector<double>(1, 1.0); }
	virtual std::vector<std::string> getSecondStageRowNames(int scen) { return std::vector<std::string>(ncons2_, "theta"); }
	virtual double scenarioProbability(int scen) { return 1.0/nscen_; }
	virtual bool isSecondStageColInteger(int scen, int col) { return false; } /** why do we need this? */

	// returns the column-oriented first-stage constraint matrix (A matrix) 
	virtual CoinPackedMatrix getFirstStageConstraints();
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen);
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen);

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
	virtual CoinPackedMatrix getFirstStageHessian();
    virtual CoinPackedMatrix getSecondStageCrossHessian(int scen);
    virtual CoinPackedMatrix getSecondStageHessian(int scen);

	/* TODO
        virtual int nLinkCons() { return 0; } 
        virtual int nLinkECons() { return 0; }
        virtual int nLinkICons() { return 0; }
        virtual CoinPackedMatrix getLinkMatrix(int nodeid){ return CoinPackedMatrix(); }
	std::string datarootname;
	int useInputDate;
	*/
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWBUNDLEPIPSINPUT_H_ */
