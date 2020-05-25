/*
 * BdSub.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "Utility/DspMessage.h"
#include "Solver/Benders/BdSub.h"
#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"

#ifdef DSP_HAS_CPX
#define DSP_OSI DspOsiCpx
#else
#define DSP_OSI DspOsiClp
#endif

BdSub::BdSub(DspParams* par):
par_(par),
nsubprobs_(0), 
subindices_(NULL),
mat_mp_(NULL), 
cglp_(NULL),
objvals_(NULL), 
solutions_(NULL),
warm_start_(NULL), 
status_(NULL) {}

BdSub::BdSub(const BdSub& rhs) :
par_(rhs.par_),
nsubprobs_(rhs.nsubprobs_) {
	/** copy things ... */
	setSubIndices(rhs.nsubprobs_, rhs.subindices_);
	mat_mp_     = new CoinPackedMatrix * [nsubprobs_];
	cglp_       = new DspOsi * [nsubprobs_];
	warm_start_ = new CoinWarmStart * [nsubprobs_];
	objvals_    = new double [nsubprobs_];
	solutions_  = new double * [nsubprobs_];
	status_     = new int [nsubprobs_];
	for (int i = 0; i < nsubprobs_; ++i) {
		mat_mp_[i] = new CoinPackedMatrix(*(rhs.mat_mp_[i]));
		cglp_[i] = rhs.cglp_[i]->clone();
		warm_start_[i] = rhs.warm_start_[i]->clone();
		objvals_[i] = rhs.objvals_[i];
		solutions_[i] = new double [cglp_[i]->si_->getNumCols()];
		CoinCopyN(rhs.solutions_[i], cglp_[i]->si_->getNumCols(), solutions_[i]);
		status_[i] = rhs.status_[i];
	}
}


BdSub::~BdSub()
{
	FREE_ARRAY_PTR(subindices_);
	FREE_2D_PTR(nsubprobs_, mat_mp_);
	FREE_2D_PTR(nsubprobs_, cglp_);
	FREE_ARRAY_PTR(objvals_);
	FREE_2D_PTR(nsubprobs_, solutions_);
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
				obj_reco, rlbd_reco, rubd_reco));

		/** creating solver interface */
		cglp_[i] = new DSP_OSI;
		cglp_[i]->si_->loadProblem(*mat_reco, clbd_reco, cubd_reco, obj_reco, rlbd_reco, rubd_reco);
		for (int j = 0; j < cglp_[i]->si_->getNumCols(); ++j)
		{
			if (ctype_reco[j] != 'C')
				cglp_[i]->si_->setInteger(j);
		}
		cglp_[i]->setLogLevel(0);

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
	int status;
	double ** Tx = NULL;

	BGN_TRY_CATCH

	/** retrieve the number of rows in subproblem */
	int nrows = mat_mp_[0]->getNumRows();

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
	FREE_PTR(cglp)

	DspOsi * cglp = NULL;

	BGN_TRY_CATCH

	DSPdebugMessage("Scenario %d\n", s);

	/** local variables */
	const double * rlbd = cgl->cglp_[s]->si_->getRowLower();
	const double * rubd = cgl->cglp_[s]->si_->getRowUpper();
	const double * clbd = cgl->cglp_[s]->si_->getColLower();
	const double * cubd = cgl->cglp_[s]->si_->getColUpper();
	const double * pi   = NULL; /** dual variables */
	const double * rc   = NULL; /** reduced costs */

	/** clone CGLP */
	cglp = new DSP_OSI;
	cglp->si_->loadProblem(*(cgl->cglp_[s]->si_->getMatrixByCol()),
			cgl->cglp_[s]->si_->getColLower(), cgl->cglp_[s]->si_->getColUpper(),
			cgl->cglp_[s]->si_->getObjCoefficients(),
			cgl->cglp_[s]->si_->getRowLower(), cgl->cglp_[s]->si_->getRowUpper());
	cglp->si_->setWarmStart(cgl->warm_start_[s]);

	int nrows = cglp->si_->getNumRows();

	/** loop over CGLP rows to update row bounds */
	for (int i = nrows - 1; i >= 0; --i)
	{
		DSPdebugMessage("s %d, i %d, rlbd %e, rubd %e, Tx %e\n", s, i, rlbd[i], rubd[i], Tx[s][i]);
		if (rlbd[i] > -1.0e+20)
			cglp->si_->setRowLower(i, rlbd[i] - Tx[s][i]);
		if (rubd[i] < 1.0e+20)
			cglp->si_->setRowUpper(i, rubd[i] - Tx[s][i]);
	}

	/** solve feasibility problem */
	int nAddedCols = 0;

	/** change to and solve feasibility problem */
	solveFeasProblem(cglp, nAddedCols);

	/** get warmstart */
	if (!cgl->warm_start_[s])
	{
		cgl->warm_start_[s] = cglp->si_->getWarmStart();
		CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(cgl->warm_start_[s]);
		basis->resize(cglp->si_->getNumRows(), cglp->si_->getNumCols() - nAddedCols);
		basis = NULL;
	}

	/** infeasible? */
	if (cglp->getPrimObjValue() > 1.e-6)
	{
		/** save solution status */
		cgl->status_[s] = DSP_STAT_PRIM_INFEASIBLE;
		DSPdebugMessage("  Subproblem %d is infeasible (infeasibility %e).\n", s, cglp->getPrimObjValue());

		/** calculate feasibility cut elements */
		calculateCutElements(
				cglp->si_->getNumRows(), cglp->si_->getNumCols(), cgl->mat_mp_[s],
				rlbd, rubd, cglp->si_->getColLower(), cglp->si_->getColUpper(),
				cglp->si_->getRowPrice(), cglp->si_->getReducedCost(), cutval[s], cutrhs[s]);

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
	chgToOrgProblem(cglp, cgl->cglp_[s]->si_->getObjCoefficients(), nAddedCols);

	/** prepare for generating optimality cuts */
#ifdef DSP_DEBUG
	cglp->setLogLevel(5);
#else
	cglp->setLogLevel(0);
#endif
	cglp->si_->setHintParam(OsiDoScale);
	cglp->si_->setHintParam(OsiDoPresolveInResolve);
	cglp->si_->setHintParam(OsiDoDualInResolve);

	/** set warmstart */
	cglp->si_->setWarmStart(cgl->warm_start_[s]);

	/** solve */
	cglp->solve();

	/** update warmstart */
	FREE_PTR(cgl->warm_start_[s]);
	cgl->warm_start_[s] = cglp->si_->getWarmStart();

	DSPdebugMessage("  objective value %E\n", cglp->getPrimObjValue());

	/** solution status */
	cgl->status_[s] = cglp->status();
	DSPdebugMessage("  solution status: %d\n", cgl->status_[s]);

	if (cgl->status_[s] == DSP_STAT_OPTIMAL)
	{
		/** TODO: add parametric cuts */

		/** get objective value */
		cgl->objvals_[s] = cglp->getPrimObjValue();

		/** get primal solution */
		CoinCopyN(cglp->si_->getColSolution(), cglp->si_->getNumCols(), cgl->solutions_[s]);

		/** get dual variables and reduced costs */
		pi = cglp->si_->getRowPrice();
		rc = cglp->si_->getReducedCost();
#ifdef DSP_DEBUG2
		for (int i = 0; i < cglp->si_->getNumRows(); ++i)
		{
			printf("  rlbd[%d] = %E, rubd[%d] = %E, pi[%d] = %E\n",
					i, rlbd[i], i, rubd[i], i, pi[i]);
		}
		for (int j = 0; j < cglp->si_->getNumCols(); ++j)
		{
			printf("  y[%d] = %E, rc[%d] = %E\n", j, cglp->si_->getColSolution()[j], j, rc[j]);
		}
#endif

		/** calculate cut elements */
		calculateCutElements(cglp->si_->getNumRows(), cglp->si_->getNumCols(),
				cgl->mat_mp_[s], rlbd, rubd, clbd, cubd, pi, rc, cutval[s], cutrhs[s]);
	} else {
		cglp->si_->writeMps("subprob");
	}

	END_TRY_CATCH(FREE_MEMORY)

	FREE_MEMORY
#undef FREE_MEMORY
}

DSP_RTN_CODE BdSub::solveFeasProblem(
		DspOsi * osi,
		int & nAddedCols)
{
	BGN_TRY_CATCH

	for (int j = 0; j < osi->si_->getNumCols(); ++j)
		osi->si_->setObjCoeff(j, 0.0);
	std::vector<int> columnStarts, rows;
	std::vector<double> elements, clbd, cubd, obj;
	for (int i = 0; i < osi->si_->getNumRows(); ++i)
	{
		columnStarts.push_back(elements.size());
		rows.push_back(i);
		if (osi->si_->getRowSense()[i] == 'G')
			elements.push_back(1.0);
		else
			elements.push_back(-1.0);
		clbd.push_back(0.0);
		cubd.push_back(COIN_DBL_MAX);
		obj.push_back(1.0);
		if (osi->si_->getRowSense()[i] == 'E' ||
			osi->si_->getRowSense()[i] == 'R')
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
	osi->si_->addCols(nAddedCols, &columnStarts[0], &rows[0], &elements[0], &clbd[0], &cubd[0], &obj[0]);

#ifdef DSP_DEBUG
	osi->setLogLevel(1);
#else
	osi->setLogLevel(0);
#endif
	osi->si_->setHintParam(OsiDoScale);
	osi->si_->setHintParam(OsiDoPresolveInInitial);
	osi->si_->setHintParam(OsiDoDualInInitial);

	osi->solve();
	//si_copy->writeMps("debug");

	/** should be always optimal */
	assert(osi->si_->isProvenOptimal());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdSub::chgToOrgProblem(
		DspOsi * osi,
		const double * obj,
		int & nAddedCols)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(delind)

	int * delind = NULL;

	BGN_TRY_CATCH

	/** delete columns */
	delind = new int [nAddedCols];
	CoinIotaN(delind, nAddedCols, osi->si_->getNumCols() - nAddedCols);
	osi->si_->deleteCols(nAddedCols, delind);

	/** set original objective */
	osi->si_->setObjective(obj);

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
