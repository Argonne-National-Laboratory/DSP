/*
 * TssBd.h
 *
 *  Created on: Dec 5, 2014
 *      Author: kibaekkim
 */

#ifndef TSSBD_H_
#define TSSBD_H_

/** DSP */
#include "SolverInterface/SolverInterface.h"
#include "Solver/TssBdSub.h"
#include "Solver/TssSolver.h"

struct TssBdStatistics
{
	double cut_generation_time_cpu_phase1_;  /**< cpu time in cut generation */
	double cut_generation_time_cpu_phase2_;  /**< cpu time in cut generation */
	double cut_generation_time_wall_phase1_; /**< wallclock time in cut generation */
	double cut_generation_time_wall_phase2_; /**< wallclock time in cut generation */
};

class TssBd: public TssSolver
{
public:

	/** default constructor */
	TssBd();

	/** default destructor */
	virtual ~TssBd();

	/** solve */
	virtual STO_RTN_CODE solve();

	/** set auxiliary variable data */
	void setAuxColData(int size, double * obj, double * clbd, double * cubd);

	/** set auxiliary variable data */
	void setAugScenarios(int size, int * indices);

private:

	/** construct master problem */
	STO_RTN_CODE constructMasterProblem(TssBdSub * tssbdsub, double lowerbound);

#if 0
	/** construct subproblems */
	STO_RTN_CODE constructSubProblem(
			double * probability,
			CoinPackedMatrix **&   mat_tech,
			OsiSolverInterface **& si_reco);
#endif

	/** to find a lower bound by solving a set of group subproblems */
	STO_RTN_CODE findLowerBound(
			const double * probability,
			double & lowerbound);

	/** configure Stochastic LP */
	STO_RTN_CODE configureSLP();

	/** solve Stochastic LP */
	STO_RTN_CODE solveSLP(TssBdSub * tssbdsub);

	/** configure Stochastic MILP */
	STO_RTN_CODE configureSMILP(TssBdSub * tssbdsub);

	/** solve Stochastic MILP */
	STO_RTN_CODE solveSMILP();

public:

	TssBdStatistics stat_;

protected:

	SolverInterface * si_; /**< my solver interface */

	int      naugs_;    /**< number of scenarios augmented to master */
	int *    augs_;     /**< augmented scenario indices */
	int      naux_;     /**< number of auxiliary variables */
	double * obj_aux_;  /**< auxiliary variable objectives */
	double * clbd_aux_; /**< auxiliary variable lower bounds */
	double * cubd_aux_; /**< auxiliary variable upper bounds */
};

#endif /* TSSBD_H_ */
