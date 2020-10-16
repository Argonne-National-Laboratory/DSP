/*
 * BdSub.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
// #define DSP_PROFILE

#include "Utility/DspUtility.h"
#include "Model/TssModel.h"
#include "Solver/Benders/BdSub.h"
#include "SolverInterface/DspOsiScip.h"
#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"

BdSub::BdSub(DspParams *par) : par_(par),
							   nsubprobs_(0),
							   subindices_(NULL),
							   probability_(NULL),
							   mat_mp_(NULL),
							   cglp_(NULL),
							   warm_start_(NULL),
							   objvals_(NULL),
							   solutions_(NULL),
							   status_(NULL),
							   recourse_has_integer_(false)
{
	names_statistics_.push_back("create ipsub");
	names_statistics_.push_back("create ipsub DspOsi");
	names_statistics_.push_back("create lpsub");
	names_statistics_.push_back("solve ipsub");
	names_statistics_.push_back("solve lpsub");
	for (unsigned i = 0; i < names_statistics_.size(); ++i)
	{
		count_statistics_[names_statistics_[i]] = 0;
		time_statistics_[names_statistics_[i]] = 0.0;
	}
}

BdSub::BdSub(const BdSub &rhs) : par_(rhs.par_),
								 nsubprobs_(rhs.nsubprobs_),
								 recourse_has_integer_(rhs.recourse_has_integer_),
								 names_statistics_(rhs.names_statistics_),
								 count_statistics_(rhs.count_statistics_),
								 time_statistics_(rhs.time_statistics_)
{
	/** copy things ... */
	setSubIndices(rhs.nsubprobs_, rhs.subindices_);
	probability_ = new double [nsubprobs_];
	mat_mp_     = new CoinPackedMatrix * [nsubprobs_];
	cglp_       = new DspOsi * [nsubprobs_];
	warm_start_ = new CoinWarmStart * [nsubprobs_];
	objvals_    = new double [nsubprobs_];
	solutions_  = new double * [nsubprobs_];
	status_     = new int [nsubprobs_];
	for (int i = 0; i < nsubprobs_; ++i) {
		probability_[i] = rhs.probability_[i];
		mat_mp_[i] = new CoinPackedMatrix(*(rhs.mat_mp_[i]));
		cglp_[i] = rhs.cglp_[i]->clone();
		cglp_[i]->setLogLevel(0);
		cglp_[i]->setNumCores(par_->getIntParam("BD/SUB/THREADS"));
		warm_start_[i] = rhs.warm_start_[i]->clone();
		objvals_[i] = rhs.objvals_[i];
		solutions_[i] = new double [cglp_[i]->si_->getNumCols()];
		CoinCopyN(rhs.solutions_[i], cglp_[i]->si_->getNumCols(), solutions_[i]);
		status_[i] = rhs.status_[i];
	}
}

BdSub::~BdSub()
{
#ifdef DSP_PROFILE
	fstream fp;
	fp.open("BdSub_statistics.txt", ios::out);
	for (unsigned i = 0; i < names_statistics_.size(); ++i)
	{
		fp << names_statistics_[i] << ","
		   << time_statistics_[names_statistics_[i]] << ","
		   << count_statistics_[names_statistics_[i]] << endl;
	}
	fp.close();
#endif

	FREE_ARRAY_PTR(subindices_);
	FREE_ARRAY_PTR(probability_);
	FREE_2D_PTR(nsubprobs_, mat_mp_);
	FREE_2D_PTR(nsubprobs_, cglp_);
	FREE_ARRAY_PTR(objvals_);
	FREE_2D_ARRAY_PTR(nsubprobs_, solutions_);
	FREE_2D_PTR(nsubprobs_, warm_start_);
	FREE_ARRAY_PTR(status_);
}

DSP_RTN_CODE BdSub::setSubIndices(int size, const int* indices)
{
	BGN_TRY_CATCH

	nsubprobs_ = size;

	/** allocate memory */
	FREE_ARRAY_PTR(subindices_);
	subindices_ = new int [size];

	CoinCopyN(indices, size, subindices_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdSub::loadProblem(DecModel* model)
{
#define FREE_MEMORY            \
	FREE_PTR(mat_reco)         \
	FREE_ARRAY_PTR(clbd_reco)  \
	FREE_ARRAY_PTR(cubd_reco)  \
	FREE_ARRAY_PTR(ctype_reco) \
	FREE_ARRAY_PTR(obj_reco)   \
	FREE_ARRAY_PTR(rlbd_reco)  \
	FREE_ARRAY_PTR(rubd_reco)

	/** recourse problem data */
	CoinPackedMatrix * mat_reco = NULL;
	double * clbd_reco   = NULL;
	double * cubd_reco   = NULL;
	double * obj_reco    = NULL;
	char *   ctype_reco  = NULL;
	double * rlbd_reco   = NULL;
	double * rubd_reco   = NULL;
	assert(model);

	BGN_TRY_CATCH

	/** allocate memory */
	probability_ = new double [nsubprobs_];
	mat_mp_     = new CoinPackedMatrix * [nsubprobs_];
	cglp_       = new DspOsi * [nsubprobs_];
	warm_start_ = new CoinWarmStart * [nsubprobs_];
	objvals_    = new double [nsubprobs_];
	solutions_  = new double * [nsubprobs_];
	status_     = new int [nsubprobs_];

	for (int i = 0; i < nsubprobs_; ++i)
	{
		/** initialize memory */
		mat_mp_[i]     = NULL;
		cglp_[i]       = NULL;
		warm_start_[i] = NULL;

		/** copy recourse problem */
		DSP_RTN_CHECK_THROW(model->copyRecoProb(subindices_[i],
				mat_mp_[i], mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, rlbd_reco, rubd_reco, false));
		
		/** get probability; 1.0 for non-stochastic */
		probability_[i] = model->isStochastic() ? dynamic_cast<TssModel*>(model)->getProbability()[i] : 1.0;

		/** creating solver interface */
		cglp_[i] = createDspOsi(par_->getIntParam("BD/SUB/SOLVER"));
		if (!cglp_[i]) throw CoinError("Failed to create DspOsi", "loadProblem", "BdSub");
		
		/** load problem */
		cglp_[i]->si_->loadProblem(*mat_reco, clbd_reco, cubd_reco, obj_reco, rlbd_reco, rubd_reco);
		for (int j = 0; j < cglp_[i]->si_->getNumCols(); ++j)
		{
			if (ctype_reco[j] != 'C') {
				cglp_[i]->si_->setInteger(j);
				recourse_has_integer_ = true;
			}
		}
		cglp_[i]->setLogLevel(0);
		cglp_[i]->setNumCores(par_->getIntParam("BD/SUB/THREADS"));

		/** allocate memory for solution */
		solutions_[i] = new double [cglp_[i]->si_->getNumCols()];

		/** free memory */
		FREE_MEMORY
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

int BdSub::generateCuts(
		int            ncols,  /**< number of first-stage variables */
		const double * x,      /**< first-stage solution */
		double **      cutval, /** dense cut coefficients for each subproblem */
		double *       cutrhs  /** cut rhs for each subproblem */)
{
#define FREE_MEMORY \
	FREE_2D_ARRAY_PTR(nsubprobs_,Tx);

	assert(mat_mp_);
	int status = DSP_STAT_UNKNOWN;
	double ** Tx = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	Tx = new double * [nsubprobs_];
	for (int s = nsubprobs_ - 1; s >= 0; --s)
	{
		Tx[s] = new double [mat_mp_[s]->getNumRows()];
		/** calculate Tx */
		mat_mp_[s]->times(x, Tx[s]);
		/** cut placeholder */
		cutval[s] = new double [ncols];
		CoinZeroN(cutval[s], ncols);
	}
	CoinZeroN(cutrhs, nsubprobs_);

	/** loop over subproblems */
	bool doContinue = true;
	for (int s = nsubprobs_ - 1; s >= 0; --s)
	{
		if (doContinue == false)
		{
			status_[s] = DSP_STAT_NOT_SOLVED;
			continue;
		}
		solveOneSubproblem(this, s, x, Tx, cutval, cutrhs);
		if (status_[s] == DSP_STAT_PRIM_INFEASIBLE)
			doContinue = false;
	}

	END_TRY_CATCH(FREE_MEMORY)

	/** free memory */
	FREE_MEMORY

	return status;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdSub::evaluateRecourse(
		const double * x, /**< [in] first-stage solution */
		double * objvals  /**< [out] objective values */) {
	double* Tx = NULL;

	for (int s = 0; s < nsubprobs_; ++s) {
		Tx = new double [mat_mp_[s]->getNumRows()];
		// calculate Tx
		mat_mp_[s]->times(x, Tx);
		// evaluate intege recourse
		DSP_RTN_CHECK_RTN_CODE(solveOneIntegerSubproblem(this, s, x, Tx, objvals));
		// free Tx
		FREE_ARRAY_PTR(Tx);
	}
	return DSP_RTN_OK;
}

void BdSub::solveOneSubproblem(
		BdSub *     cgl,
		int            s,            /**< scenario index */
		const double * x,            /**< first-stage solution */
		double **      Tx,           /**< Tx */
		double **      cutval,       /**< Benders cut body */
		double *       cutrhs,       /**< Benders cut RHS */
		int            enableOptCuts /**< whether to generate optimality cuts or not */)
{
#define FREE_MEMORY \
	FREE_PTR(osi)

	DspOsi *osi = NULL;
	OsiSolverInterface *si = NULL;

	BGN_TRY_CATCH

	DSPdebugMessage("Scenario %d\n", s);

	/** local variables */
	const double * rlbd = cgl->cglp_[s]->si_->getRowLower();
	const double * rubd = cgl->cglp_[s]->si_->getRowUpper();
	const double * clbd = cgl->cglp_[s]->si_->getColLower();
	const double * cubd = cgl->cglp_[s]->si_->getColUpper();
	const double * pi   = NULL; /** dual variables */
	const double * rc   = NULL; /** reduced costs */
	const char * ctype = cgl->cglp_[s]->si_->getColType();

	double stime = CoinGetTimeOfDay(); // tic

	/** clone CGLP */
	osi = cgl->cglp_[s]->clone();
	osi->setLogLevel(0);
	osi->setNumCores(cgl->par_->getIntParam("BD/SUB/THREADS"));
	si = osi->si_;
	si->setWarmStart(cgl->warm_start_[s]);

	int nrows = si->getNumRows();

	/** loop over CGLP rows to update row bounds */
	for (int i = nrows - 1; i >= 0; --i)
	{
		DSPdebugMessage("s %d, i %d, rlbd %e, rubd %e, Tx %e\n", s, i, rlbd[i], rubd[i], Tx[s][i]);
		if (rlbd[i] > -1.0e+20)
			si->setRowLower(i, rlbd[i] - Tx[s][i]);
		if (rubd[i] < 1.0e+20)
			si->setRowUpper(i, rubd[i] - Tx[s][i]);
	}

	/** collect statistics */
	cgl->count_statistics_["create lpsub"]++;
	cgl->time_statistics_["create lpsub"] += CoinGetTimeOfDay() - stime; // toc

	/** solve feasibility problem */
	int nAddedCols = 0;

	/** change to and solve feasibility problem */
	solveFeasProblem(si, nAddedCols);

	/** get warmstart */
	if (!cgl->warm_start_[s])
	{
		cgl->warm_start_[s] = si->getWarmStart();
		if (cgl->warm_start_[s] != NULL)
		{
			CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(cgl->warm_start_[s]);
			basis->resize(si->getNumRows(), si->getNumCols() - nAddedCols);
			basis = NULL;
		}
	}

	/** infeasible? */
	if (si->getObjValue() > 1.e-6)
	{
		/** save solution status */
		cgl->status_[s] = DSP_STAT_PRIM_INFEASIBLE;
		DSPdebugMessage("  Subproblem %d is infeasible (infeasibility %e).\n", s, si->getObjValue());

		/** calculate feasibility cut elements */
		calculateCutElements(
			si->getNumRows(), si->getNumCols(), cgl->mat_mp_[s],
			rlbd, rubd, si->getColLower(), si->getColUpper(),
			si->getRowPrice(), si->getReducedCost(), cutval[s], cutrhs[s]);

		FREE_MEMORY
		return;
	}

	if (!enableOptCuts)
	{
		cgl->status_[s] = DSP_STAT_FEASIBLE;
		FREE_MEMORY
		return;
	}
	/** change to original problem */
	chgToOrgProblem(si, cgl->cglp_[s]->si_->getObjCoefficients(), nAddedCols);

	/** prepare for generating optimality cuts */
#ifdef DSP_DEBUG
	si->messageHandler()->setLogLevel(5);
#else
	si->messageHandler()->setLogLevel(0);
#endif
	si->setHintParam(OsiDoScale);
	si->setHintParam(OsiDoPresolveInResolve);
	si->setHintParam(OsiDoDualInResolve);

	/** set warmstart */
	si->setWarmStart(cgl->warm_start_[s]);

	stime = CoinGetTimeOfDay(); // tic

	/** solve */
	si->initialSolve();

	/** collect statistics */
	cgl->count_statistics_["solve lpsub"]++;
	cgl->time_statistics_["solve lpsub"] += CoinGetTimeOfDay() - stime; // toc

	/** update warmstart */
	FREE_PTR(cgl->warm_start_[s]);
	cgl->warm_start_[s] = si->getWarmStart();

	DSPdebugMessage("  objective value %E\n", si->getObjValue());

	/** solution status */
	cgl->status_[s] = DspOsi::dsp_status(si);
	DSPdebugMessage("  solution status: %d\n", cgl->status_[s]);

	if (cgl->status_[s] == DSP_STAT_OPTIMAL)
	{
		/** TODO: add parametric cuts */

		/** get objective value */
		cgl->objvals_[s] = si->getObjValue();

		/** get primal solution */
		CoinCopyN(si->getColSolution(), si->getNumCols(), cgl->solutions_[s]);

		/** get dual variables and reduced costs */
		pi = si->getRowPrice();
		rc = si->getReducedCost();
#ifdef DSP_DEBUG2
		for (int i = 0; i < si->getNumRows(); ++i)
		{
			printf("  rlbd[%d] = %E, rubd[%d] = %E, pi[%d] = %E\n",
					i, rlbd[i], i, rubd[i], i, pi[i]);
		}
		for (int j = 0; j < si->getNumCols(); ++j)
		{
			printf("  y[%d] = %E, rc[%d] = %E\n", j, si->getColSolution()[j], j, rc[j]);
		}
#endif

		/** calculate cut elements */
		calculateCutElements(si->getNumRows(), si->getNumCols(),
							 cgl->mat_mp_[s], rlbd, rubd, clbd, cubd, pi, rc, cutval[s], cutrhs[s]);
	} else {
		si->writeMps("subprob");
	}

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY
#undef FREE_MEMORY
}

DSP_RTN_CODE BdSub::solveOneIntegerSubproblem(
			BdSub *     cgl,
			int            s,     /**< scenario index */
			const double * x,     /**< first-stage solution */
			double *       Tx,    /**< Tx */
			double *       objval /**< objective value */)
{
#define FREE_MEMORY \
	FREE_PTR(osi)

	DspOsi *osi = NULL;
	OsiSolverInterface *si = NULL;
	DSP_RTN_CODE ret = DSP_RTN_OK;

	BGN_TRY_CATCH

	DSPdebugMessage("Scenario %d\n", s);

	/** local variables */
	const double * rlbd = cgl->cglp_[s]->si_->getRowLower();
	const double * rubd = cgl->cglp_[s]->si_->getRowUpper();
	const char * ctype = cgl->cglp_[s]->si_->getColType();

	double stime = CoinGetTimeOfDay(); // tic

	osi = cgl->cglp_[s]->clone();
	osi->setLogLevel(0);
	osi->setNumCores(cgl->par_->getIntParam("BD/SUB/THREADS"));
	si = osi->si_;

	/** collect statistics */
	cgl->count_statistics_["create ipsub DspOsi"]++;
	cgl->time_statistics_["create ipsub DspOsi"] += CoinGetTimeOfDay() - stime; // toc
	stime = CoinGetTimeOfDay();

	/** mark integer variables */
	for (int j = 0; j < si->getNumCols(); ++j)
	{
		if (ctype[j] != 'C')
			si->setInteger(j);
	}

	/** loop over CGLP rows to update row bounds */
	for (int i = 0; i < si->getNumRows(); ++i)
	{
		if (rlbd[i] > -1.0e+20)
			si->setRowLower(i, rlbd[i] - Tx[i]);
		if (rubd[i] < 1.0e+20)
			si->setRowUpper(i, rubd[i] - Tx[i]);
	}

	/** quite */
	si->messageHandler()->setLogLevel(0);

	/** collect statistics */
	cgl->count_statistics_["create ipsub"]++;
	cgl->time_statistics_["create ipsub"] += CoinGetTimeOfDay() - stime; // toc
	stime = CoinGetTimeOfDay();											 // tic

	/** solve */
	si->initialSolve();
	if (si->getNumIntegers() > 0)
		si->branchAndBound();
	DSPdebugMessage("  objective value %E\n", si->getObjValue());

	/** collect statistics */
	cgl->count_statistics_["solve ipsub"]++;
	cgl->time_statistics_["solve ipsub"] += CoinGetTimeOfDay() - stime; // toc

	/** solution status */
	cgl->status_[s] = DspOsi::dsp_status(si);
	DSPdebugMessage("  solution status: %d\n", cgl->status_[s]);

	if (cgl->status_[s] == DSP_STAT_OPTIMAL) {
		/** get objective value */
		objval[s] = si->getObjValue();
	} else {
		printf("Unexpected solution status: s %d status %d\n", s, cgl->status_[s]);
		objval[s] = 1.0e+20;
		si->writeMps("int_subprob");
		ret = DSP_RTN_ERR;
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return ret;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdSub::solveFeasProblem(
	OsiSolverInterface *si,
	int &nAddedCols)
{
	BGN_TRY_CATCH

	for (int j = 0; j < si->getNumCols(); ++j)
		si->setObjCoeff(j, 0.0);
	std::vector<int> columnStarts, rows;
	std::vector<double> elements, clbd, cubd, obj;
	for (int i = 0; i < si->getNumRows(); ++i)
	{
		columnStarts.push_back(elements.size());
		rows.push_back(i);
		if (si->getRowSense()[i] == 'G')
			elements.push_back(1.0);
		else
			elements.push_back(-1.0);
		clbd.push_back(0.0);
		cubd.push_back(COIN_DBL_MAX);
		obj.push_back(1.0);
		if (si->getRowSense()[i] == 'E' ||
			si->getRowSense()[i] == 'R')
		{
			columnStarts.push_back(elements.size());
			rows.push_back(i);
			elements.push_back(1.0);
			clbd.push_back(0.0);
			cubd.push_back(COIN_DBL_MAX);
			obj.push_back(1.0);
		}
	}
	columnStarts.push_back(elements.size());
	nAddedCols = columnStarts.size() - 1;
	si->addCols(nAddedCols, &columnStarts[0], &rows[0], &elements[0], &clbd[0], &cubd[0], &obj[0]);

#ifdef DSP_DEBUG
	si->messageHandler()->setLogLevel(1);
#else
	si->messageHandler()->setLogLevel(0);
#endif
	si->setHintParam(OsiDoScale);
	si->setHintParam(OsiDoPresolveInInitial);
	si->setHintParam(OsiDoDualInInitial);

	si->initialSolve();
	//si_copy->writeMps("debug");

	/** should be always optimal */
	assert(si->isProvenOptimal());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdSub::chgToOrgProblem(
	OsiSolverInterface *si,
	const double *obj,
	int &nAddedCols)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(delind)

	int * delind = NULL;

	BGN_TRY_CATCH

	/** delete columns */
	delind = new int [nAddedCols];
	CoinIotaN(delind, nAddedCols, si->getNumCols() - nAddedCols);
	si->deleteCols(nAddedCols, delind);

	/** set original objective */
	si->setObjective(obj);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE BdSub::calculateCutElements(
		int nrows,                         /**< [in] number of rows in subproblem */
		int ncols,                         /**< [in] number of columns in subproblem */
		const CoinPackedMatrix * mat_tech, /**< [in] technology matrix */
		const double * rlbd,               /**< [in] row lower bounds */
		const double * rubd,               /**< [in] row upper bounds */
		const double * clbd,               /**< [in] column lower bounds */
		const double * cubd,               /**< [in] column upper bounds */
		const double * pi,                 /**< [in] dual variables corresponding to constraints */
		const double * rc,                 /**< [in] reduced cost corresponding to variables */
		double *       cutval,             /**< [out] cut coefficients */
		double &       cutrhs              /**< [out] cut rhs */)
{
	BGN_TRY_CATCH

	int i, j;

	/** calculate pi^T T */
	mat_tech->transposeTimes(pi, cutval);

	/** calculate pi^T h */
	/** loop over rows */
	for (i = nrows - 1; i >= 0; --i)
	{
		//if (fabs(pi[i]) < 1.0e-10) continue;
        if (rubd[i] >= COIN_DBL_MAX)
        {
            /** Ax >= b */
        	cutrhs += pi[i] * rlbd[i];
        }
        else if (rlbd[i] <= -COIN_DBL_MAX)
        {
            /** Ax <= b */
        	cutrhs += pi[i] * rubd[i];
        }
        else if (rlbd[i] == rubd[i])
        {
            /** Ax = b */
        	cutrhs += pi[i] * rlbd[i];
        }
        else if (pi[i] >= 0.0)
        {
            /** l <= Ax <= u, bounded by l (-y + s_L - s_U = 0) */
        	cutrhs += pi[i] * rlbd[i];
        }
        else
        {
            /** l <= Ax <= u, bounded by u */
        	cutrhs += pi[i] * rubd[i];
        }
		//DSPdebugMessage("cutrhs %e, pi %e rlbd %e rubd %e\n", cutrhs, pi[i], rlbd[i], rubd[i]);
	}
	//DSPdebugMessage("cutrhs %E\n", cutrhs);

	/** loop over columns */
	for (j = ncols - 1; j >= 0; --j)
	{
		//if (fabs(rc[j]) < 1.0e-10) continue;
		//cutval[j] += rc[j];
        if (cubd[j] >= COIN_DBL_MAX)
        {
        	/** free variable? */
        	if (clbd[j] <= -COIN_DBL_MAX)
        		continue;
            /** x_j >= l_j */
            cutrhs += rc[j] * clbd[j];
        }
        else if (clbd[j] <= -COIN_DBL_MAX)
        {
            /** x_j <= u_j */
            cutrhs += rc[j] * cubd[j];
        }
        else if (clbd[j] == cubd[j])
        {
            /** l_j = x_j = u_j */
            cutrhs += rc[j] * clbd[j];
        }
        else if (rc[j] >= 0.0)
        {
            /** l_j <= x_j <= u_j, bounded by l */
            cutrhs += rc[j] * clbd[j];
        }
        else
        {
            /** l_j <= x_j <= u_j, bounded by u */
            cutrhs += rc[j] * cubd[j];
        }
		//DSPdebugMessage("cutrhs %e, rc %e clbd %e cubd %e\n", cutrhs, rc[j], clbd[j], cubd[j]);
	}

	//DSPdebugMessage("cutrhs %E\n", cutrhs);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DspOsi * BdSub::createDspOsi(int solver) {
	DspOsi * osi = NULL;
	BGN_TRY_CATCH

	switch (solver) {
		case OsiScip:
			osi = new DspOsiScip();
			break;
		case OsiCpx:
#ifdef DSP_HAS_CPX
			osi = new DspOsiCpx();
#else
			throw CoinError("Cplex is not available.", "createDspOsi", "BdSub");
#endif
			break;
		case OsiGrb:
#ifdef DSP_HAS_GRB
			osi = new DspOsiGrb();
#else
			throw CoinError("Gurobi is not available.", "createDspOsi", "BdSub");
#endif
			break;
		case OsiClp:
			osi = new DspOsiClp();
			break;
		default:
			char coinmsg[128];
			sprintf(coinmsg, "Invalid parameter value (solver = %d)", solver);
			throw CoinError(coinmsg, "createDspOsi", "BdSub");
			break;
	}

	END_TRY_CATCH_RTN(;,osi);
	return osi;
}
