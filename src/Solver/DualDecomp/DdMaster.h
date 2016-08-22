/*
 * DdMaster.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMASTER_H_
#define SRC_SOLVER_DUALDECOMP_DDMASTER_H_

#include "Solver/DualDecomp/DdSolver.h"
#include "SolverInterface/SolverInterface.h"

class DdMaster: public DdSolver {

	friend class DdMW;
	friend class DdMWSerial;
	friend class DdMWAsync;
	friend class DdMWSync;

public:

	/** constructor */
	DdMaster(
			DspParams *  par,    /**< parameter pointer */
			DecModel *   model,  /**< model pointer */
			DspMessage * message /**< message pointer */);
			//int nworkers          /**< number of workers */);

	/** destructor */
	virtual ~DdMaster();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** update problem */
	virtual DSP_RTN_CODE updateProblem() = 0;

	/** set init solution */
	virtual DSP_RTN_CODE setInitSolution(const double * sol);

	/** termination test */
	virtual DSP_RTN_CODE terminationTest() {return status_;}

public:

	const double * getLambda() {return lambda_;}

	SolverInterface * getSiPtr() {return si_;}

	double   getBestPrimalObjective() {return bestprimobj_;}
	double   getBestDualObjective()   {return bestdualobj_;}
	double * getBestPrimalSolution()  {return bestprimsol_;}
	double * getBestDualSolution()    {return bestdualsol_;}
	double getAbsDualityGap() {return fabs(bestprimobj_-bestdualobj_);}
	double getRelDualityGap() {return fabs(bestprimobj_-bestdualobj_) / (1.0e-10 + fabs(bestprimobj_));}
	double getAbsApproxGap() {return fabs(primobj_-bestdualobj_);}
	double getRelApproxGap() {return fabs(primobj_-bestdualobj_) / (1.0e-10 + fabs(primobj_));}

protected:

	SolverInterface * si_; /**< solver interface */

	double bestprimobj_;   /**< best primal objective value */
	double bestdualobj_;   /**< best dual objective value */
	double * bestprimsol_; /** best primal solution */
	double * bestdualsol_; /** best dual solution */
	double * lambda_; /**< lambda part of the solution */
//	int nworkers_;       /**< number of workers */

	//int nsubprobs_;         /**< number of subproblems for the current worker */
	//int * subindex_;        /**< array of subproblem indices */
	double * subprimobj_;   /**< subproblem primal objective values */
	double * subdualobj_;   /**< subproblem dual objective values */
	double ** subsolution_; /**< subproblem solution */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMASTER_H_ */
