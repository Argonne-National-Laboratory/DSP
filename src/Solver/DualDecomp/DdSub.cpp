/*
 * DdSub.h
 *
 *  Renamed on: Feb 16, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

/** Coin */
#include "OsiCbcSolverInterface.hpp"
#ifdef DSP_HAS_CPLEX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
#endif
/** Dsp */
#include "Utility/DspMacros.h"
#include "Utility/DspMessage.h"
#include "Solver/DualDecomp/DdSub.h"
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#include "SolverInterface/SCIPbranchruleLB.h"

/** default constructor */
DdSub::DdSub(int s, DspParams * par, DecModel * model, DspMessage * message) :
	DecSolver(model, par, message),
	si_(NULL),
	sind_(s),
	nrows_coupling_(0),
	ncols_coupling_(0),
	theta_(-COIN_DBL_MAX),
	gapTol_(0.0001),
	obj_(NULL),
	lambda_(NULL),
	cpl_mat_(NULL),
	cpl_cols_(NULL),
	cpl_rhs_(NULL),
	obj_offset_(0),
	parRelaxIntegrality_(NULL)
{
	/** nothing to do */
}

DdSub::~DdSub()
{
	FREE_PTR(si_);
	FREE_PTR(obj_);
	FREE_ARRAY_PTR(lambda_);
	FREE_ARRAY_PTR(cpl_rhs_);
	obj_offset_ = 0;
	parRelaxIntegrality_ = NULL;
}

DSP_RTN_CODE DdSub::init()
{
	BGN_TRY_CATCH

	/** create problem */
	createProblem();

	/** allocate memory */
	primsol_ = new double [si_->getNumCols()];

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** solve problem */
DSP_RTN_CODE DdSub::solve()
{
	si_->resolve();

	/** check status. there might be unexpected results. */
	if (si_->isProvenOptimal()) {
		primobj_ = si_->getObjValue();
		/** TODO: how to get dual bound? */
		dualobj_ = si_->getObjValue();
		CoinCopyN(si_->getColSolution(), si_->getNumCols(), primsol_);
		DSPdebugMessage("primal objective %+e\n", primobj_);
		status_ = DSP_STAT_OPTIMAL;
	} else if (si_->isIterationLimitReached())
		status_ = DSP_STAT_LIM_ITERorTIME;
	else if (si_->isAbandoned())
		status_ = DSP_STAT_ABORT;
	else if (si_->isProvenPrimalInfeasible())
		status_ = DSP_STAT_PRIM_INFEASIBLE;
	else if (si_->isProvenDualInfeasible())
		status_ = DSP_STAT_DUAL_INFEASIBLE;

	return DSP_RTN_OK;
}

/** create problem */
DSP_RTN_CODE DdSub::createProblem() {
#define FREE_MEMORY       \
    FREE_PTR(mat)         \
    FREE_ARRAY_PTR(clbd)  \
    FREE_ARRAY_PTR(cubd)  \
    FREE_ARRAY_PTR(ctype) \
    FREE_ARRAY_PTR(rlbd)  \
    FREE_ARRAY_PTR(rubd)

    /** problem data */
    CoinPackedMatrix *mat = NULL;
    double *clbd = NULL;
    double *cubd = NULL;
    char *ctype = NULL;
    double *rlbd = NULL;
    double *rubd = NULL;
    int augs[1];
    double clbd_aux[1];
    double cubd_aux[1];
    double obj_aux[1];
    int cpl_ncols;

    BGN_TRY_CATCH

        /** parameters */
        parRelaxIntegrality_ = par_->getBoolPtrParam("RELAX_INTEGRALITY");
        gapTol_ = par_->getDblParam("SCIP/GAP_TOL");

        /** augmented subproblem index */
        augs[0] = sind_;

        /** for auxiliary term */
        clbd_aux[0] = -COIN_DBL_MAX;
        cubd_aux[0] = COIN_DBL_MAX;
        obj_aux[0] = 0.0;

        /** decompose model */
        DSP_RTN_CHECK_THROW(
                model_->decompose(1, augs, 1, clbd_aux, cubd_aux, obj_aux,
                                  mat, clbd, cubd, ctype, obj_, rlbd, rubd));

        DSP_RTN_CHECK_THROW(
                model_->decomposeCoupling(1, augs, cpl_mat_, cpl_cols_, cpl_ncols));

        /** number of coupling variables and constraints for this subproblem */
        ncols_coupling_ = model_->getNumSubproblemCouplingCols(sind_);
        nrows_coupling_ = model_->getNumSubproblemCouplingRows(sind_);
        assert(ncols_coupling_ == cpl_ncols);
        assert(nrows_coupling_ == cpl_mat_->getNumRows());

        /** storage for lambda */
        lambda_ = new double[nrows_coupling_];
        CoinZeroN(lambda_, nrows_coupling_);

        /** copy right-hand side of coupling rows */
        cpl_rhs_ = new double[nrows_coupling_];
        for (int i = 0; i < nrows_coupling_; i++)
            cpl_rhs_[i] = model_->getRhsCouplingRow(i);
        obj_offset_ = 0;

        /** adjust first-stage cost */
        if (model_->isStochastic()) {
            TssModel *tssModel;
            try {
                tssModel = dynamic_cast<TssModel *>(model_);
            }
            catch (const std::bad_cast &e) {
                printf("Error: Model claims to be stochastic when it is not\n");
                return DSP_RTN_ERR;
            }

            double probability = tssModel->getProbability()[sind_];
            for (int j = 0; j < tssModel->getNumCols(0); ++j)
                obj_[j] *= probability;

            /** convert column types */
            if (parRelaxIntegrality_[0]) {
                for (int j = 0; j < tssModel->getNumCols(0); ++j)
                    ctype[j] = 'C';
            }
            if (parRelaxIntegrality_[1])
                CoinFillN(ctype + tssModel->getNumCols(0), tssModel->getNumCols(1), 'C');
        } else {
            if (parRelaxIntegrality_[0] || parRelaxIntegrality_[1]) {
                for (int j = 0; j < mat->getNumCols(); j++)
                    ctype[j] = 'C';
            }
        }

#ifdef DSP_HAS_CPLEX
        si_ = new OsiCpxSolverInterface();
#else
        si_ = new OsiCbcSolverInterface();
#endif

        /** no display */
        si_->messageHandler()->logLevel(0);

        /** load problem */
        si_->loadProblem(*mat, clbd, cubd, obj_, rlbd, rubd);
        DSPdebug(mat->verifyMtx(4));

        /** set integers */
        for (int j = 0; j < si_->getNumCols(); ++j) {
        	if (ctype[j] != 'C')
        		si_->setInteger(j);
        }

        /** set solution gap tolerance */
        setGapTol(gapTol_);

    END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

    FREE_MEMORY

    return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** update problem */
DSP_RTN_CODE DdSub::updateProblem(
		double * lambda,
		double primal_bound)
{
#define FREE_MEMORY       \
	FREE_ARRAY_PTR(newobj);

	double * newobj = NULL;

	BGN_TRY_CATCH

	int ncols = si_->getNumCols();

	if (lambda)
	{
		/** copy lambda */
		CoinCopyN(lambda, nrows_coupling_, lambda_);

		/** allocate memory */
		newobj = new double [ncols];
		obj_offset_ = 0;

		/** update objective coefficients */
		assert(obj_);
		for (int j = 0; j < ncols; j++)
			newobj[j] = obj_[j];
		for (int i = 0; i < nrows_coupling_; i++)
		{
			/** add lambda wrt coupling row i */
			assert(!cpl_mat_->isColOrdered()); /** matrix must be by row */
			int size = cpl_mat_->getVector(i).getNumElements();
			const int * inds = cpl_mat_->getVector(i).getIndices();
			const double * elems = cpl_mat_->getVector(i).getElements();
			for (int j = 0; j < size; j++)
				newobj[inds[j]] += lambda[i] * elems[j];

			/* if rhs is not zero, then the objective has a constant offset, which is added when passing to the master */
			obj_offset_ += lambda[i] * (-cpl_rhs_[i]);
		}

		si_->setObjective(newobj);
	}

	/** update primal bound (bounds of auxiliary constraint) */
	if (primal_bound < COIN_DBL_MAX)
		si_->setColUpper(ncols - 1, primal_bound);//si_->setColBounds(ncols - 1, primal_bound, primal_bound);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** push cuts */
DSP_RTN_CODE DdSub::pushCuts(OsiCuts * cuts)
{
	si_->applyCuts(*cuts);
	return DSP_RTN_OK;
}

/** set wall clock time limit */
void DdSub::setTimeLimit(double sec)
{
	assert(si_);
#ifdef DSP_HAS_CPLEX
	OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
	if (cpx) CPXsetdblparam(cpx->getEnvironmentPtr(), CPX_PARAM_TILIM, sec);
#else
	OsiCbcSolverInterface* cbc = dynamic_cast<OsiCbcSolverInterface*>(si_);
	if (cbc) cbc->setMaximumSeconds(sec);
#endif
}

/** set accuracy tolerance */
void DdSub::setGapTol(double tol)
{
	assert(si_);
	gapTol_ = tol;
#ifdef DSP_HAS_CPLEX
	OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
	if (cpx) CPXsetdblparam(cpx->getEnvironmentPtr(), CPX_PARAM_EPGAP, gapTol_);
#else
	OsiCbcSolverInterface* cbc = dynamic_cast<OsiCbcSolverInterface*>(si_);
	if (cbc) cbc->getModelPtr()->setAllowableFractionGap(gapTol_);
#endif
}

/** set print level */
void DdSub::setPrintLevel(int level) {
	si_->messageHandler()->logLevel(level);
}
