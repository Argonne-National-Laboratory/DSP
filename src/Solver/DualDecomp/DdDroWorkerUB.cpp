// #define DSP_DEBUG
#include "Model/TssModel.h"
#include "Solver/DualDecomp/DdDroWorkerUB.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"

#ifdef DSP_HAS_SCIP
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#endif

DdDroWorkerUB::DdDroWorkerUB(const DdDroWorkerUB &rhs) : DdWorkerUB(rhs)
{
	osi_dro_ = rhs.osi_dro_->clone();
}

DdDroWorkerUB::~DdDroWorkerUB()
{
	FREE_PTR(osi_dro_);
}

DSP_RTN_CODE DdDroWorkerUB::init()
{
	BGN_TRY_CATCH
	status_ = DSP_STAT_MW_CONTINUE;

	/** Create the stochastic upper bounding subproblems.
	 * We still need to solve these subproblems,
	 * in addition to the DRO upper bounding problem.
	 */
	DSP_RTN_CHECK_THROW(DdWorkerUB::createProblem());

	/** Create the DRO upper bounding problem. */
	DSP_RTN_CHECK_THROW(createProblem());

	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DdDroWorkerUB::createProblem()
{
#define FREE_MEMORY          \
	FREE_ARRAY_PTR(clbd_dro) \
	FREE_ARRAY_PTR(cubd_dro) \
	FREE_ARRAY_PTR(rlbd_dro) \
	FREE_ARRAY_PTR(rubd_dro) \
	FREE_ARRAY_PTR(obj_dro)  \
	FREE_ARRAY_PTR(elem_dro) \
	FREE_ARRAY_PTR(ind_dro)  \
	FREE_ARRAY_PTR(bgn_dro)  \
	FREE_ARRAY_PTR(len_dro)

	/** DRO UB problem */
	double *clbd_dro = NULL;
	double *cubd_dro = NULL;
	double *obj_dro = NULL;
	double *rlbd_dro = NULL;
	double *rubd_dro = NULL;
	double *elem_dro = NULL;
	int *ind_dro = NULL;
	int *bgn_dro = NULL;
	int *len_dro = NULL;

	TssModel *tss = NULL;

	BGN_TRY_CATCH

	try
	{
		tss = dynamic_cast<TssModel *>(model_);
	}
	catch (const std::bad_cast &e)
	{
		printf("This is not a stochastic programming problem.\n");
		return DSP_RTN_ERR;
	}

	/** create DRO upper bounding problem */
	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");
	int ncols_dro = 1 + model_->getNumReferences();
	int nrows_dro = model_->getNumReferences() * tss->getNumScenarios();
	int nzcnt_dro = nrows_dro * 2;

	/** The second-stage objective coefficients of each subproblem 
	 * need re-sacled by the original probability. 
	 */
	for (int s = 0; s < nsubprobs; ++s)
	{
		const double *obj_reco = osi_[s]->si_->getObjCoefficients();
		for (int j = 0; j < tss->getNumCols(1); ++j)
		{
			osi_[s]->si_->setObjCoeff(j, obj_reco[j] / tss->getProbability()[par_->getIntPtrParam("ARR_PROC_IDX")[s]]);
		}
	}

	/** allocate memory */
	clbd_dro = new double[ncols_dro];
	cubd_dro = new double[ncols_dro];
	obj_dro = new double[ncols_dro];
	rlbd_dro = new double[nrows_dro];
	rubd_dro = new double[nrows_dro];
	bgn_dro = new int[nrows_dro + 1];
	len_dro = new int[nrows_dro];
	ind_dro = new int[nzcnt_dro];
	elem_dro = new double[nzcnt_dro];

	/** column bounds */
	clbd_dro[0] = 0.0;
	CoinFillN(clbd_dro + 1, ncols_dro - 1, -COIN_DBL_MAX);
	CoinFillN(cubd_dro, ncols_dro, COIN_DBL_MAX);

	/** row bounds */
	CoinFillN(rlbd_dro, nrows_dro, -COIN_DBL_MAX);
	CoinFillN(rubd_dro, nrows_dro, COIN_DBL_MAX);

	/** objective coefficients */
	obj_dro[0] = model_->getWassersteinSize();
	for (int j = 1; j < ncols_dro; ++j)
		obj_dro[j] = model_->getReferenceProbability(j - 1);

	int pos_dro = 0;
	int rnum = 0;
	for (int k = 0; k < tss->getNumScenarios(); ++k)
	{
		for (int s = 0; s < model_->getNumReferences(); ++s)
		{
			rnum = k * model_->getNumReferences() + s;

			bgn_dro[rnum] = pos_dro;

			// alpha
			ind_dro[pos_dro] = 0;
			elem_dro[pos_dro] = model_->getWassersteinDist(s, k);
			pos_dro++;

			// beta_s
			ind_dro[pos_dro] = 1 + s;
			elem_dro[pos_dro] = 1.0;
			pos_dro++;

			len_dro[rnum] = pos_dro - bgn_dro[rnum];
		}
	}
	bgn_dro[nrows_dro] = pos_dro;
	assert(pos_dro == nzcnt_dro);

	/** constraint matrix */
	CoinPackedMatrix *mat_dro = new CoinPackedMatrix(false, ncols_dro, nrows_dro, nzcnt_dro, elem_dro, ind_dro, bgn_dro, len_dro);
	DSPdebug(mat_dro->verifyMtx(4));

	/** creating solver interface */
	osi_dro_ = createDspOsi();
	if (!osi_dro_)
		throw CoinError("Failed to create DspOsi", "createProblem", "DdWorkerUB");

	/** no display */
	osi_dro_->setLogLevel(0);

	/** load problem */
	osi_dro_->si_->loadProblem(*mat_dro, clbd_dro, cubd_dro, obj_dro, rlbd_dro, rubd_dro);

	/** free matrix */
	FREE_PTR(mat_dro);

	END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdDroWorkerUB::solve()
{
	double cputime;
	double walltime;

	BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime = 0.0;
	double total_walltime = 0.0;
	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	for (unsigned s = 0; s < nsubprobs; ++s)
	{
		cputime = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		/** set time limit */
		osi_[s]->setTimeLimit(
			CoinMin(CoinMax(0.01, time_remains_),
					par_->getDblParam("DD/SUB/TIME_LIM")));

		/** solve */
		osi_[s]->solve();

		/** check status. there might be unexpected results. */
		int status = osi_[s]->status();
		DSPdebugMessage("status = %d\n", status);
		switch (status)
		{
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
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
			primobj = COIN_DBL_MAX;
			dualobj = -COIN_DBL_MAX;
			break;
		}

		primobj += osi_[s]->getPrimObjValue();
		dualobj += osi_[s]->getDualObjValue();
		CoinCopyN(osi_[s]->si_->getColSolution(), osi_[s]->si_->getNumCols(), &primsols_[s][0]);
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** consume time */
		time_remains_ -= CoinGetTimeOfDay() - walltime;
	}

	/** solve the DRO UB problem */
	if (model_->isDro() && status_ != DSP_STAT_MW_STOP)
	{
		assert(nsubprobs == model_->getNumSubproblems());
		for (int k = 0; k < nsubprobs; ++k)
		{
			for (int s = 0; s < model_->getNumReferences(); ++s)
			{
				osi_dro_->si_->setRowLower(k * model_->getNumReferences() + s, osi_[k]->si_->getObjValue());
			}
		}
		osi_dro_->solve();

		if (osi_dro_->si_->isProvenOptimal())
		{
			primobj = osi_dro_->getPrimObjValue();
		}
		else
		{
			primobj = COIN_DBL_MAX;
		}
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
