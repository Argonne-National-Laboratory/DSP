/*
 * DdSub.h
 *
 *  Renamed on: Feb 16, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include <Utility/DspMacros.h>
#include <Utility/DspMessage.h>
#include "Model/TssModel.h"
#include "Solver/DualDecomp/DdSub.h"
#include "SolverInterface/SolverInterfaceClp.h"

#ifndef NO_CPX
	#include "SolverInterface/SolverInterfaceCpx.h"
#endif

#ifndef NO_SCIP
	#include "SolverInterface/SolverInterfaceScip.h"
	#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
	#include "SolverInterface/SCIPbranchruleLB.h"
#endif

/** default constructor */
DdSub::DdSub(int s, DspParams * par, DecModel * model, DspMessage * message) :
	DecSolver(par, model, message),
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
	DSP_RTN_CHECK_THROW(createProblem());

	/** add cut generator (lazycuts) */
	DSP_RTN_CHECK_THROW(addCutGenerator());

	/** allocate memory */
	primsol_ = new double [si_->getNumCols()];

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** solve problem */
DSP_RTN_CODE DdSub::solve()
{
	si_->solve();

	/** check status. there might be unexpected results. */
	status_ = si_->getStatus();
	switch (status_) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_LIM_ITERorTIME:
	case DSP_STAT_STOPPED_GAP:
	case DSP_STAT_STOPPED_NODE:
	case DSP_STAT_STOPPED_TIME:
		primobj_ = si_->getPrimalBound();
		dualobj_ = si_->getDualBound();
		CoinCopyN(si_->getSolution(), si_->getNumCols(), primsol_);
		DSPdebugMessage("primal objective %+e\n", primobj_);
		break;
	default:
		break;
	}

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
    gapTol_ = par_->getDblParam("MIP/GAP_TOL");

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

    /** number of integer variables */
    int nIntegers = model_->getNumIntegers();

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
            for (int j = 0; j < tssModel->getNumCols(0); ++j) {
                if (ctype[j] != 'C')
                    nIntegers--;
                ctype[j] = 'C';
            }
        }
        if (parRelaxIntegrality_[1]) {
            for (int j = 0; j < tssModel->getNumCols(1); ++j) {
                if (ctype[tssModel->getNumCols(0) + j] != 'C')
                    nIntegers--;
            }
            CoinFillN(ctype + tssModel->getNumCols(0), tssModel->getNumCols(1), 'C');
        }
    } else {
        if (parRelaxIntegrality_[0] || parRelaxIntegrality_[1]) {
            for (int j = 0; j < mat->getNumCols(); j++) {
                if (ctype[j] != 'C')
                    nIntegers--;
                ctype[j] = 'C';
            }
        }
    }

    if (nIntegers > 0) {
    	switch (par_->getIntParam("SOLVER/MIP")) {
    	case CPLEX:
#ifndef NO_CPX
    		si_ = new SolverInterfaceCpx(par_);
    		break;
#endif
    	case EXT_SCIP:
#ifndef NO_SCIP
            si_ = new SolverInterfaceScip(par_);
            break;
#endif
    	default:
    		break;
    	}
    } else
        si_ = new SolverInterfaceClp(par_);

    /** no display */
    si_->setPrintLevel(0);

    /** load problem */
    si_->loadProblem(mat, clbd, cubd, obj_, ctype, rlbd, rubd, "DdSub");
    DSPdebug(mat->verifyMtx(4));

    /** set solution gap tolerance */
    si_->setGapTol(gapTol_);

    END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

    FREE_MEMORY

    return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** add cut generator */
DSP_RTN_CODE DdSub::addCutGenerator() {
#ifndef NO_SCIP
    SolverInterfaceScip *si = dynamic_cast<SolverInterfaceScip *>(si_);
    if (si) {
        /** create constraint handler */
        SCIPconshdlrBendersDd *conshdlr = new SCIPconshdlrBendersDd(si->getSCIP());
        conshdlr->setOriginalVariables(si->getNumCols(), si->getSCIPvars(), 1);
        DSPdebugMessage("numcols %d, conshdlr %p\n", si->getNumCols(), conshdlr);
        /** add constraint handler */
        si->addConstraintHandler(conshdlr, true, true);
    } else {
        /** TODO */
        printf("Warning: Cut generation supports only SCIP.\n");
    }
#endif
    return DSP_RTN_OK;
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

		si_->setObjCoef(newobj);
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
#ifndef NO_SCIP
	/** TODO only for SCIP */
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		si->setCuts(cuts);
	}
#endif
	return DSP_RTN_OK;
}

/** set wall clock time limit */
void DdSub::setTimeLimit(double sec)
{
	assert(si_);
	si_->setTimeLimit(sec);
}

/** set accuracy tolerance */
void DdSub::setGapTol(double tol)
{
	assert(si_);
	gapTol_ = tol;
	si_->setGapTol(gapTol_);
}

/** set print level */
void DdSub::setPrintLevel(int level)
{
#ifndef NO_SCIP
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
		si_->setPrintLevel(level);
	else
		si_->setPrintLevel(CoinMax(0, level - 2));
#endif
}
