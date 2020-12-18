/**
 * DspOsiCpx.h
 *
 * 12/12/2019
 * Kibaek Kim
 */
  
#ifndef SRC_SOLVERINTERFACE_DSPOSICPX_H_
#define SRC_SOLVERINTERFACE_DSPOSICPX_H_

// #define DSP_DEBUG

#ifdef DSP_HAS_CPX

#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#include "SolverInterface/DspOsi.h"
#include "Utility/DspMessage.h"

class DspOsiCpx : public DspOsi {
public:

	/** default constructor */
	DspOsiCpx() : 
	DspOsi()
	{
		si_ = new OsiCpxSolverInterface();
		cpx_ = dynamic_cast<OsiCpxSolverInterface*>(si_);
	}

	/** copy constructor */
	DspOsiCpx(const DspOsiCpx& rhs) : 
	DspOsi(rhs) {
        si_ = new OsiCpxSolverInterface(*(rhs.cpx_));
		cpx_ = dynamic_cast<OsiCpxSolverInterface*>(si_);
	}

	/** clone constructor */
	virtual DspOsiCpx* clone() const {
		return new DspOsiCpx(*this);
	}

	/** destructor */
	virtual ~DspOsiCpx() {
		cpx_ = NULL;
	}

	/** load quadratic objective */
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {
		if (mat.isColOrdered()) {
			for (int j = 0; j < mat.getMajorDim(); ++j) {
				for (int k = 0; k < mat.getVectorSize(j); ++k) {
					int i = mat.getIndices()[mat.getVectorStarts()[j] + k];
					double v = mat.getElements()[mat.getVectorStarts()[j] + k];
					if (i==j){
						CPXchgqpcoef(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
						i, j, 2*v);
					}
					else{
						CPXchgqpcoef(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
						i, j, v);
					}
					
				}
			}
		} else {
			for (int i = 0; i < mat.getMajorDim(); ++i) {
				for (int k = 0; k < mat.getVectorSize(i); ++k) {
					int j = mat.getIndices()[mat.getVectorStarts()[i] + k];
					double v = mat.getElements()[mat.getVectorStarts()[i] + k];
					if (i==j){
						CPXchgqpcoef(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
						i, j, 2*v);
					}
					else{
						CPXchgqpcoef(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL),
						i, j, v);
					}
				}
			}
		}
		isqp_ = true;
	}

	
	/* For a continuous problem:
	 * (1) If there is a quadratic term in the objective function and all the constraints in the model are linear, the problem is termed a QP.
	 * * Default: the barrier optimizer.
	 * * can be changed to primal/dual/network simplex.
	 * * (For more information about optimizing a QP, see the website: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/UsrMan/topics/cont_optim/qp/11_optimize.html) 
	 * (2) If the model has any constraints containing a quadratic term, regardless of the objective function, the problem is termed a QCP.
	 * * The barrier optimizer is the only optimizer available to solve QCPs.
	 * * (For more information about solving a QCP, see the website: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/UsrMan/topics/cont_optim/qcp/15_solve.html) 
	 * For a discrete problem:
	 * (1) If there is a quadratic term in the objective function and all the constraints in the model are linear, the problem is termed a MIQP.
	 * * (For more information about solving a MIQP, see the website: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/UsrMan/topics/discr_optim/mip_quadratic/02_introMIQP.html) 
	 * (2) If the model has any constraints containing a quadratic term, regardless of the objective function, the problem is termed a MIQCP
	 * * (For more information about solving a MIQCP, see the website: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/UsrMan/topics/discr_optim/mip_quadratic/03_introMIQCP.html) 
	 */
	/** load quadratic constrs */
	virtual void addQuadraticRows(int nqrows, int * linnzcnt, int * quadnzcnt, double * rhs, int * sense, int ** linind, double ** linval, int ** quadrow, int ** quadcol, double ** quadval)
	{
		if (nqrows > 0)
			isqcp_ = true;
		for (int i = 0; i < nqrows; i++) 
		{
			int err = CPXaddqconstr(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), linnzcnt[i], quadnzcnt[i], rhs[i], sense[i], 
				linind[i], linval[i], quadrow[i], quadcol[i], quadval[i], NULL);
			checkDspOsiError(err, "CPXaddqconstr", "addQuadraticRows");	
		}
	}

	/** throw error */
	virtual inline void checkDspOsiError(int err, std::string cpxfuncname, std::string dsposimethod) 
	{
		if (err != 0) {
			char s[100];
			sprintf(s, "%s returned error %d", cpxfuncname.c_str(), err);
			throw CoinError(s, dsposimethod.c_str(), "DspOsiCpx");
		}
	}

	/** write problem file */
	virtual void writeProb(char const * filename_str, char const * filetype_str)
	{
		int err = CPXwriteprob(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), filename_str, filetype_str);
		checkDspOsiError(err, "CPXwriteprob", "writeProb");
	}

	/** change problem type to MIQCP */
	virtual void switchToMIQCP(void) {
		DSPdebugMessage("DspOsiCpx::switchToMIQCP()\n");	

  		CPXENVptr env = cpx_->getEnvironmentPtr();
		CPXLPptr lp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

		int err = CPXchgprobtype(env, lp, CPXPROB_MIQCP);
    	checkDspOsiError(err, "CPXchgprobtype", "switchToMIQCP");

		int nc = cpx_->getNumCols();
    	int *cindarray = new int[nc];

    	for (int i = 0; i < nc; ++i)
      		cindarray[i] = i;

		err = CPXchgctype(env, lp, nc, cindarray, cpx_->getCtype());
    	checkDspOsiError(err, "CPXchgctype", "switchToMIQCP");

    	delete[] cindarray;
	}

	/** change problem type to MIQP */
	virtual void switchToMIQP(void) {

		DSPdebugMessage("DspOsiCpx::switchToMIQP()\n");

  		CPXENVptr env = cpx_->getEnvironmentPtr();
		CPXLPptr lp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

		int err = CPXchgprobtype(env, lp, CPXPROB_MIQP);
    	checkDspOsiError(err, "CPXchgprobtype", "switchToMIQP");

    	int nc = cpx_->getNumCols();
    	int *cindarray = new int[nc];

    	for (int i = 0; i < nc; ++i)
      		cindarray[i] = i;

		err = CPXchgctype(env, lp, nc, cindarray, cpx_->getCtype());
    	checkDspOsiError(err, "CPXchgctype", "switchToMIQP");

    	delete[] cindarray;
  	}

	/** change problem type to QCP */
	virtual void switchToQCP(void)
	{
		DSPdebugMessage("DspOsiCpx::switchToQCP()\n");

		CPXENVptr env = cpx_->getEnvironmentPtr();
		CPXLPptr lp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

		int err = CPXchgprobtype(env, lp, CPXPROB_QCP);
		checkDspOsiError(err, "CPXchgprobtype", "switchToQCP");

	}

	/** change problem type to QP */
	virtual void switchToQP(void)
	{
		DSPdebugMessage("DspOsiCpx::switchToQP()\n");

  		CPXENVptr env = cpx_->getEnvironmentPtr();
		CPXLPptr lp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

		int err = CPXchgprobtype(env, lp, CPXPROB_QP);
		checkDspOsiError(err, "CPXchgprobtype", "switchToQP");
  		
	}

	/** solve problem */
	virtual void solve() 
	{
		if (!isqp_ && !isqcp_) {
			if (si_->getNumIntegers() > 0) {
				DSPdebugMessage("DspOsiCpx::solve(), si_->branchAndBound() \n");
				ismip_ = true;
				si_->branchAndBound();
			} else {
				DSPdebugMessage("DspOsiCpx::solve(), si_->initialSolve() \n");
				ismip_ = false;
				si_->initialSolve();
			}
		}
		else 
		{
			if (si_->getNumIntegers() > 0) {
				DSPdebugMessage("DspOsiCpx::solve(), MIQCQPSolve() \n");
				ismip_ = true;
				MIQCQPSolve();
			} else {
				DSPdebugMessage("DspOsiCpx::solve(), QCQPSolve() \n");
				ismip_ = false;
				QCQPSolve();
			}	
		}
	}

	/** modified OsiCpxSolverInterface::initialSolve */
	void QCQPSolve()
	{
  		if (isqp_)
		  	switchToQP();
		else if (isqcp_)
			switchToQCP();

  		bool takeHint;
  		OsiHintStrength strength;

  		int algorithm = 0;
		/* Only barrier method is available for solving QCP */
  		if (!isqcp_) 
		{
			cpx_->getHintParam(OsiDoDualInInitial, takeHint, strength);
  			if (strength != OsiHintIgnore)
    			algorithm = takeHint ? -1 : 1;
		}

  		int presolve = 1;
  		cpx_->getHintParam(OsiDoPresolveInInitial, takeHint, strength);
  		if (strength != OsiHintIgnore)
    		presolve = takeHint ? 1 : 0;

		CPXENVptr env = cpx_->getEnvironmentPtr();
  		CPXLPptr lp = cpx_->getLpPtr(OsiCpxSolverInterface::FREECACHED_RESULTS);

  		if (presolve)
    		CPXsetintparam(env, CPX_PARAM_PREIND, CPX_ON);
  		else
    		CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);

		if (cpx_->messageHandler()->logLevel() == 0)
    		CPXsetintparam(env, CPX_PARAM_SIMDISPLAY, 0);
  		else if (cpx_->messageHandler()->logLevel() == 1)
    		CPXsetintparam(env, CPX_PARAM_SIMDISPLAY, 1);
  		else if (cpx_->messageHandler()->logLevel() > 1)
    		CPXsetintparam(env, CPX_PARAM_SIMDISPLAY, 2);

		/* Only barrier method is available for solving QCP */
		if (!isqcp_) 
		{	
			/** need access to the private member disableadvbasis of OsiCpxSolverInterface 
			 * CPXsetintparam(env, CPX_PARAM_ADVIND, !cpx_->disableadvbasis);
			 * Just keep default setting of CPX_PARAM_ADVIND
			 */
  			
			double objoffset;
			double primalobjlimit;
			double dualobjlimit;
			cpx_->getDblParam(OsiObjOffset, objoffset);
			cpx_->getDblParam(OsiPrimalObjectiveLimit, primalobjlimit);
			cpx_->getDblParam(OsiDualObjectiveLimit, dualobjlimit);

			if (cpx_->getObjSense() == +1) {
				if (primalobjlimit < COIN_DBL_MAX)
				CPXsetdblparam(env, CPX_PARAM_OBJLLIM, primalobjlimit + objoffset);
				if (dualobjlimit > -COIN_DBL_MAX)
				CPXsetdblparam(env, CPX_PARAM_OBJULIM, dualobjlimit + objoffset);
			} else {
				if (primalobjlimit > -COIN_DBL_MAX)
				CPXsetdblparam(env, CPX_PARAM_OBJULIM, primalobjlimit + objoffset);
				if (dualobjlimit < COIN_DBL_MAX)
				CPXsetdblparam(env, CPX_PARAM_OBJLLIM, dualobjlimit + objoffset);
			}
		}

		int term;
		bool resolved = false;
		while (1) {
			switch (algorithm) {
			default:
			case 0:
				if (isqcp_) 
				{
					term = CPXbaropt(env, lp);
					checkDspOsiError(term, "CPXqcpopt", "QCQPSolve");
				}
				else if (isqp_)
				{
					term = CPXqpopt(env, lp);
					checkDspOsiError(term, "CPXqpopt", "QCQPSolve");
				}
				break;
			case 1:
				if (isqcp_) 
				{
					throw CoinError("Primal simplex is not available for QCP.", "QCQPSolve", "DspOsiCpx");
				}
				else if (isqp_)
				{
					CPXsetintparam(env, CPX_PARAM_QPMETHOD, 1);
					term = CPXqpopt(env, lp);
					checkDspOsiError(term, "CPXqpopt", "QCQPSolve");
				}
				break;
			case -1:
				if (isqcp_) 
				{
					throw CoinError("Dual simplex is not available for QCP.", "QCQPSolve", "DspOsiCpx");
				}
				else if (isqp_)
				{
					CPXsetintparam(env, CPX_PARAM_QPMETHOD, 2);
					term = CPXqpopt(env, lp);
					checkDspOsiError(term, "CPXqpopt", "QCQPSolve");
				}
				break;
			}

			/* If the problem is found infeasible during presolve, resolve it to get a proper term code */
			if (resolved) {
				CPXsetintparam(env, CPX_PARAM_PREIND, CPX_ON);
				break;
			}
			else 
			{
				int stat = CPXgetstat(env, cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
				if (stat == CPX_STAT_INForUNBD && presolve) 
				{
					CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
				} else 
					break;
				resolved = true;
			}
		}
	}
	/** modified OsiCpxSolverInterface::branchAndBound */
	void MIQCQPSolve()
	{
  		int term;

  		if (isqp_)
		  	switchToMIQP();
		else if (isqcp_)
			switchToMIQCP();
#ifdef DSP_DEBUG
		char filename[128];
		sprintf(filename, "miqcqp_model"); 
		writeProb(filename, "lp");
#endif
		CPXENVptr env = cpx_->getEnvironmentPtr();

		/** disabled domipstart in OsiCpxSolverInterface and have the default setting of CPX_PARAM_ADVIND
		 * that was able to a column solution to CPLEX before starting MIP solve (copymipstart)
		 * in order to enable it, we need to have access to private members of OsiCpxSolverInterface (ncols_ and domipstart)
		 */

  		CPXLPptr lp = cpx_->getLpPtr(OsiCpxSolverInterface::FREECACHED_RESULTS);

		if (cpx_->messageHandler()->logLevel() == 0)
			CPXsetintparam(env, CPX_PARAM_SIMDISPLAY, 0);
		else if (cpx_->messageHandler()->logLevel() == 1)
			CPXsetintparam(env, CPX_PARAM_SIMDISPLAY, 1);
		else if (cpx_->messageHandler()->logLevel() > 1)
			CPXsetintparam(env, CPX_PARAM_SIMDISPLAY, 2);

  		term = CPXmipopt(env, lp);
  		checkDspOsiError(term, "CPXmipopt", "MIQCQPSolve");
	}
	

	virtual void use_simplex() {
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
	}

	virtual void use_barrier() {
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		// CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_BARCROSSALG, -1); // This has been deprecated.
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);
		// CPXsetdblparam(cpx_->getEnvironmentPtr(), CPXPARAM_Barrier_ConvergeTol, 1e-5);
		// CPXsetintparam(cpx_->getEnvironmentPtr(), CPXPARAM_Emphasis_Numerical, CPX_ON);
		// CPXsetintparam(cpx_->getEnvironmentPtr(), CPXPARAM_Preprocessing_Reduce, 1);
		// CPXsetintparam(cpx_->getEnvironmentPtr(), CPXPARAM_Preprocessing_Presolve, 1);
	}

	/** solution statue */
	virtual int status() {
		int status = DSP_STAT_UNKNOWN;
		int probtype = CPXgetprobtype(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
		int stat = CPXgetstat(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));

		if (!ismip_) {
			switch(stat) {
			case CPX_STAT_OPTIMAL:
			case CPX_STAT_OPTIMAL_INFEAS:
				status = DSP_STAT_OPTIMAL;
				break;
			case CPX_STAT_INFEASIBLE:
				status = DSP_STAT_PRIM_INFEASIBLE;
				break;
			case CPX_STAT_UNBOUNDED:
				status = DSP_STAT_DUAL_INFEASIBLE;
				break;
			case CPX_STAT_ABORT_OBJ_LIM:
			case CPX_STAT_ABORT_PRIM_OBJ_LIM:
				status = DSP_STAT_LIM_PRIM_OBJ;
				break;
			case CPX_STAT_ABORT_DUAL_OBJ_LIM:
				status = DSP_STAT_LIM_DUAL_OBJ;
				break;
			case CPX_STAT_ABORT_IT_LIM:
				status = DSP_STAT_STOPPED_ITER;
				break;
			case CPX_STAT_ABORT_TIME_LIM:
				status = DSP_STAT_STOPPED_TIME;
				break;
			case CPX_STAT_NUM_BEST:
			case CPX_STAT_FEASIBLE:
				status = DSP_STAT_FEASIBLE;
				break;
			case CPX_STAT_ABORT_USER:
				status = DSP_STAT_STOPPED_USER;
				break;
			default:
				status = DSP_STAT_UNKNOWN;
				break;
			}
		} else if (ismip_) {
			switch(stat) {
			case CPXMIP_OPTIMAL:
			case CPXMIP_OPTIMAL_TOL:
			case CPXMIP_OPTIMAL_INFEAS:
				status = DSP_STAT_OPTIMAL;
				break;
			case CPXMIP_INFEASIBLE:
			case CPXMIP_NODE_LIM_INFEAS:
			case CPXMIP_TIME_LIM_INFEAS:
				status = DSP_STAT_PRIM_INFEASIBLE;
				break;
			case CPXMIP_UNBOUNDED:
			case CPXMIP_INForUNBD:
				status = DSP_STAT_DUAL_INFEASIBLE;
				break;
			case CPXMIP_NODE_LIM_FEAS:
				status = DSP_STAT_STOPPED_NODE;
				break;
			case CPXMIP_TIME_LIM_FEAS:
				status = DSP_STAT_STOPPED_TIME;
				break;
			case CPXMIP_FAIL_FEAS:
			case CPXMIP_FAIL_INFEAS:
			case CPXMIP_MEM_LIM_FEAS:
			case CPXMIP_MEM_LIM_INFEAS:
			case CPXMIP_ABORT_FEAS:
			case CPXMIP_ABORT_INFEAS:
			case CPXMIP_FAIL_FEAS_NO_TREE:
			case CPXMIP_FAIL_INFEAS_NO_TREE:
				status = DSP_STAT_ABORT;
				break;
			default:
				status = DSP_STAT_UNKNOWN;
				break;
			}
		}
		return status;
	}

	/** get primal objective value */
	virtual double getPrimObjValue() {
		double objval = 0.0;
		int err;
  		int solntype;
		
		CPXsolninfo(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), NULL, &solntype, NULL, NULL);
    
		if (solntype != CPX_NO_SOLN) {
      		err = CPXgetobjval(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &objval);
      		checkDspOsiError(err, "CPXgetobjval", "getObjValue");
		} else {
			// return 0.0 as objective value if no availible solution
			objval = 0.0;
		}
		return objval;
	}

	/** get dual objective value */
	virtual double getDualObjValue() {
		double val;
		if (si_->getNumIntegers() > 0) {
			CPXgetbestobjval(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &val);
			val *= cpx_->getObjSense();
		} else {
			val = si_->getObjValue();
		}
		return val;
	}

	/** get number of branch-and-bound nodes explored */
	virtual int getNumNodes() {
		if (si_->getNumIntegers() > 0)
			return CPXgetnodecnt(cpx_->getEnvironmentPtr(), cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
		else
			return 0;
	}

	/** set number of cores */
	virtual void setNumCores(int num) {
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_THREADS, CoinMax(0, num));
	}

	/** set time limit */
	virtual void setTimeLimit(double time) {
		CPXsetdblparam(cpx_->getEnvironmentPtr(), CPX_PARAM_TILIM, CoinMax(CoinMin(1.0e+75, time), 1.0e-4));
	}

	/** set node limit */
	virtual void setNodeLimit(int num)
	{
		CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_NODELIM, CoinMax(1, CoinMin(2100000000, num)));
	}

	/** set relative MIP gap */
	virtual void setRelMipGap(double tol) {
		CPXsetdblparam(cpx_->getEnvironmentPtr(), CPX_PARAM_EPGAP, CoinMax(0.0, CoinMin(1.0, tol)));
	}

    OsiCpxSolverInterface* cpx_;   
};

#endif

#endif /* SRC_SOLVERINTERFACE_DSPOSICPX_H_ */