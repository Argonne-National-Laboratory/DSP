/*
 * DdWorkerUB.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
// #define DSP_DEBUG_WRITE
#include "Model/DecTssModel.h"
#include "Solver/DualDecomp/DdWorkerUB2.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"

#ifdef DSP_HAS_SCIP
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#endif

DdWorkerUB2::DdWorkerUB2(
	DecModel *model, /**< model pointer */
	DspParams *par,	 /**< parameter pointer */
	DspMessage *message /**< message pointer */) : DdWorkerUB(model, par, message)
	// ,
												//    bestub_(COIN_DBL_MAX),
												//    ub_(0.0)
{
}

DdWorkerUB2::DdWorkerUB2(const DdWorkerUB2 &rhs) : DdWorkerUB(rhs)
// ,
												// bestub_(rhs.bestub_),
												// primsols_(rhs.primsols_),
												// ub_(rhs.ub_)
{}

DdWorkerUB2::~DdWorkerUB2()
{
	// check whether subprobs_[s] is being deleted
	// FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), mat_mp_);
	// FREE_2D_ARRAY_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), obj_org_);
	// FREE_2D_ARRAY_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rlbd_org_);
	// FREE_2D_ARRAY_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rubd_org_);
	// FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), osi_);
}

DSP_RTN_CODE DdWorkerUB2::init()
{
	DSPdebugMessage("initiating DdWorkerUB2\n");
	BGN_TRY_CATCH
	/** status */
	status_ = DSP_STAT_MW_CONTINUE;
	/** create problem */
	DSP_RTN_CHECK_THROW(createProblem(
		par_->getIntPtrParamSize("ARR_PROC_IDX"), par_->getIntPtrParam("ARR_PROC_IDX")));
	
	primsols_.resize(subprobs_.size());
	
	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		/** allocate array size for each scenario primal solution */
		primsols_[s].resize(subprobs_[s]->getDspOsiPtr()->si_->getNumCols());

		// cout << "num rows: " << subprobs_[s]->getDspOsiPtr()->si_->getNumRows() << endl;

		/** add auxiliary constraints for coupling
		 * the following auxiliary constraints are added at the end of the upper bound problem:
		 * 	x_1 = 0
		 * 	...
		 * 	x_n = 0
		 */
		int n = model_->getNumSubproblemCouplingCols(s);
		// cout << "num added rows: " << n << endl;
		int ccnt = 0; /* An integer that specifies the number of new columns in the constraints being added to the constraint matrix */
		int nrows = n;
		int nznt = n; /* An integer that specifies the number of nonzero constraint coefficients to be added to the constraint matrix. This specifies the length of the arrays rmatind and rmatval. */
		char * sense = new char[n];
		int * rmatbeg = new int[n];
		int * rmatind = new int[n];/* index for coupling vars */
		double * rmatval = new double[n];
		double * rhs = new double[n];
		for (int i = 0; i < n; i++)
		{
			rmatind[i] = model_->getSubproblemCouplingColIndices(s)[i]; 
			sense[i] = 'E';
			rmatbeg[i] = i;
			rmatval[i] = 1; 
			rhs[i] = 0; 
		}
		subprobs_[s]->getDspOsiPtr()->addRows(ccnt, nrows, nznt, rhs, sense, rmatbeg, rmatind, rmatval);
		// cout << "num updated rows: " << subprobs_[s]->getDspOsiPtr()->si_->getNumRows() << endl;

		// delete arrays
	}

	#ifdef DSP_DEBUG_WRITE
		for (unsigned s = 0; s < subprobs_.size(); ++s)
		{
			/* write in lp file */
			char filename[128];
			sprintf(filename, "DdWorkerUB2_scen%d", s); 
			DSPdebugMessage("writing initial upper bound subproblem for scenario %d in filename.lp\n", s, filename);
			subprobs_[s]->getDspOsiPtr()->writeProb(filename, "lp");
		}
	#endif
	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DdWorkerUB2::createProblem(int nsubprobs, int* subindex)
{
BGN_TRY_CATCH

	/** release all the subproblem data */
	for (unsigned s = 0; s < subprobs_.size(); ++s)
		FREE_PTR(subprobs_[s]);
	subprobs_.clear();

	for (int s = 0; s < nsubprobs; ++s) {

        /** create subproblem instance */
		/* relax integrality for first-stage variables, since they will be fixed as some scalars when evaluated */
		par_->setBoolPtrParam("RELAX_INTEGRALITY", 0, true);
        DdSub *subprob = new DdSub(subindex[s], par_, model_, message_);
		
        /** initialize */
        DSP_RTN_CHECK_THROW(subprob->init());
        assert(subprob->getSiPtr());
        /** store */
        subprobs_.push_back(subprob);
    }
	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}


double DdWorkerUB2::evaluate(int n, double *solution)
{
	std::vector<int> indices;
	std::vector<double> elements;
	for (int i = 0; i < n; ++i)
		if (fabs(solution[i]) > 1.0e-10)
		{
			indices.push_back(i);
			elements.push_back(solution[i]);
		}

	CoinPackedVector *s = new CoinPackedVector(indices.size(), &indices[0], &elements[0]);
	double ub = evaluate(s);

	delete s;

	return ub;
}

double DdWorkerUB2::evaluate(CoinPackedVector *solution)
{
#define FREE_MEMORY    \
	FREE_ARRAY_PTR(indices) \
	FREE_ARRAY_PTR(values)

	TssModel *tss = NULL;

	// BGN_TRY_CATCH

	tss = dynamic_cast<TssModel *>(model_);
	if (tss == NULL)
		throw "This is not a stochastic programming problem.";

	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		/** add auxiliary constraints for coupling */
		// int cnt = model_->getNumSubproblemCouplingCols(s);
		int n = model_->getNumSubproblemCouplingCols(s);
		int cnt = solution->getNumElements();
		/* An integer that specifies the total number of righthand side coefficients to be changed, and thus specifies the length of the arrays indices and values.*/
		int * indices = new int[cnt];
		/* An array of length cnt containing the numeric indices of the rows corresponding to the linear constraints for which righthand side coefficients are to be changed. */
		double * values = new double[cnt];
		/* An array of length cnt containing the new values of the righthand side coefficients of the linear constraints present in indices.*/
		for (int i = 0; i < cnt; i++)
		{
			/* variable indice that will be fixed */
			int var_i = solution->getIndices()[i];
			DSPdebugMessage("update upper bound subproblem %d: fix x%d at %f\n", s, var_i, solution->getElements()[i]);
			/* row indice for x_var_i = hat x_var_i */
			indices[i] = subprobs_[s]->getDspOsiPtr()->si_->getNumRows() - (n-var_i); 
			/* hat x_var_i */
			values[i] = solution->getElements()[i];
		}
		subprobs_[s]->getDspOsiPtr()->chgRhs(cnt, indices, values);
	
		#ifdef DSP_DEBUG_WRITE
				/* write in lp file */
				char filename[128];
				sprintf(filename, "DdWorkerUB_scen_updated%d", s); 
				DSPdebugMessage("writing updated upper bound subproblem for scenario %d in %s.lp\n", s, filename);
				subprobs_[s]->getDspOsiPtr()->writeProb(filename, "lp");
		#endif

		FREE_MEMORY
		// END_TRY_CATCH_RTN(FREE_MEMORY, COIN_DBL_MAX)
		
	}

	// delete arrays

	/** set initial status */
	status_ = DSP_STAT_MW_CONTINUE;

	/** solve upper bounding problem */
	DSP_RTN_CHECK_RTN_CODE(solve());

	/** update upper bound */
	if (status_ != DSP_STAT_MW_CONTINUE)
		ub_ = COIN_DBL_MAX;

	// double ub = evaluate(s);

	

	return ub_;

#undef FREE_MEMORY
	
}


DSP_RTN_CODE DdWorkerUB2::solve()
{
	double cputime;
	double walltime;

	BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime = 0.0;
	double total_walltime = 0.0;

	for (unsigned s = 0; s < subprobs_.size(); ++s)
	{
		cputime = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		/** set time limit */
		subprobs_[s]->getDspOsiPtr()->setTimeLimit(
			CoinMin(CoinMax(0.01, time_remains_),
					par_->getDblParam("DD/SUB/TIME_LIM")));

		/** solve */
		subprobs_[s]->getDspOsiPtr()->solve();

		/** check status. there might be unexpected results. */
		int status = subprobs_[s]->getDspOsiPtr()->status();
		DSPdebugMessage("status = %d\n", status);
		switch (status)
		{
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
		case DSP_STAT_FEASIBLE:
			break;
		default:
			status_ = DSP_STAT_MW_STOP;
			message_->print(10,
							"Warning: subproblem %d solution status is %d\n", s,
							status);
			break;
		}
		if (status_ == DSP_STAT_MW_STOP)
		{
			DSPdebugMessage("status_ (dsp) = DSP_STAT_MW_STOP_, status (osi) =%d\n", status);
			primobj = COIN_DBL_MAX;
			dualobj = -COIN_DBL_MAX;
			break;
		}

		primobj += subprobs_[s]->getDspOsiPtr()->getPrimObjValue();
		dualobj += subprobs_[s]->getDspOsiPtr()->getDualObjValue();
		CoinCopyN(subprobs_[s]->getDspOsiPtr()->si_->getColSolution(), subprobs_[s]->getDspOsiPtr()->si_->getNumCols(), &primsols_[s][0]);
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** consume time */
		time_remains_ -= CoinGetTimeOfDay() - walltime;
	}

	/** get primal objective */
	ub_ = primobj;
	DSPdebugMessage("ub_ = %e\n", ub_);
	DSPdebugMessage("status_ %d\n", status_);

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)

	return DSP_RTN_OK;
}
