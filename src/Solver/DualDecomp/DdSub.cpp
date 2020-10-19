/*
 * DdSub.h
 *
 *  Renamed on: Feb 16, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "Utility/DspMacros.h"
#include "Utility/DspMessage.h"
#include "Model/TssModel.h"
#include "Solver/DualDecomp/DdSub.h"
#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"

#ifdef DSP_HAS_SCIP
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#endif

/** default constructor */
DdSub::DdSub(int s, DspParams * par, DecModel * model, DspMessage * message) :
DecSolver(model, par, message),
sind_(s),
nrows_coupling_(0),
ncols_coupling_(0),
theta_(-COIN_DBL_MAX),
gapTol_(0.0001),
obj_(NULL),
cpl_mat_(NULL),
cpl_cols_(NULL),
cpl_rhs_(NULL),
obj_offset_(0),
parRelaxIntegrality_(NULL) {}

DdSub::DdSub(const DdSub& rhs) :
DecSolver(rhs),
sind_(rhs.sind_),
nrows_coupling_(rhs.nrows_coupling_),
ncols_coupling_(rhs.ncols_coupling_),
theta_(rhs.theta_),
gapTol_(rhs.gapTol_),
obj_offset_(rhs.obj_offset_) {
	// copy obj_
	obj_ = new double [getSiPtr()->getNumCols()];
	CoinCopyN(rhs.obj_, getSiPtr()->getNumCols(), obj_);

	// copy coupling matrix
	cpl_mat_ = new CoinPackedMatrix(*(rhs.cpl_mat_));

	// copy coupling columns
	cpl_cols_ = new int [ncols_coupling_];
	CoinCopyN(rhs.cpl_cols_, ncols_coupling_, cpl_cols_);

	// copy coupling rhs
	cpl_rhs_ = new double [nrows_coupling_];
	CoinCopyN(rhs.cpl_rhs_, nrows_coupling_, cpl_rhs_);
}

DdSub::~DdSub() {
	FREE_ARRAY_PTR(obj_);
	FREE_PTR(cpl_mat_);
	FREE_ARRAY_PTR(cpl_cols_);
	FREE_ARRAY_PTR(cpl_rhs_);
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
    primsol_.resize(getSiPtr()->getNumCols());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** solve problem */
DSP_RTN_CODE DdSub::solve()
{
	/** check dual infeasibility */
	bool dualinfeas = false;

	while (1) {
		osi_->solve();
	
		/** check status. there might be unexpected results. */
		status_ = osi_->status();
		DSPdebugMessage("solution status %d\n", status_);
		switch (status_) {
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
			primobj_ = osi_->getPrimObjValue();
			dualobj_ = osi_->getDualObjValue();
			assert(primsol_.size() == getSiPtr()->getNumCols());
			CoinCopyN(getSiPtr()->getColSolution(), getSiPtr()->getNumCols(), &primsol_[0]);
			DSPdebugMessage("primal objective %+e, dual objective %+e\n", primobj_, dualobj_);
			DSPdebugMessage("fraction gap %e\n", fabs(primobj_-dualobj_) / fabs(dualobj_));
			dualinfeas = false;
#ifdef DSP_DEBUG
			char submps[64];
			sprintf(submps, "sub%d", sind_);
			getSiPtr()->writeMps(submps);
#endif
			break;
		case DSP_STAT_LIM_INFEAS:
			primobj_ = COIN_DBL_MAX;
			dualobj_ = osi_->getDualObjValue();
			dualinfeas = false;
			break;
		case DSP_STAT_DUAL_INFEASIBLE:
			message_->print(0, "Subproblem %d is dual infeasible. DSP will fix any unbounded column bounds to a large number.\n", sind_);
			for (int j = 0; j < getSiPtr()->getNumCols(); ++j) {
				if (getSiPtr()->getColLower()[j] < -1.0e+20) {
					DSPdebugMessage("Fix column %d lower bound %+e to %+e\n", j, getSiPtr()->getColLower()[j], -1.0e+10);
					getSiPtr()->setColLower(j, -1.0e+10);
				}
				if (getSiPtr()->getColUpper()[j] > 1.0e+20) {
					DSPdebugMessage("Fix column %d upper bound %+e to %+e\n", j, getSiPtr()->getColUpper()[j], 1.0e+10);
					getSiPtr()->setColUpper(j, 1.0e+10);
				}
			}
			if (dualinfeas) {
				char submps[64];
				sprintf(submps, "dual_infeas_sub%d", sind_);
				getSiPtr()->writeMps(submps);
				dualinfeas = false;
			} else {
				dualinfeas = true;
			}
			//DSPdebugMessage("Dual infeasible subproblem!\n");
			//DSPdebug(getSiPtr()->writeMps("dual_infeas_sub"));
			break;
		default:
			break;
		}
		if (!dualinfeas) break;
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
    FREE_ARRAY_PTR(rubd)  \
    FREE_ARRAY_PTR(obj)

    /** problem data */
    CoinPackedMatrix *mat = NULL;
    double *clbd = NULL;
    double *cubd = NULL;
    char *ctype = NULL;
    double *rlbd = NULL;
    double *rubd = NULL;
	double *obj = NULL;
    int augs[1];
    double clbd_aux[1];
    double cubd_aux[1];
    double obj_aux[1];
    int cpl_ncols;

    BGN_TRY_CATCH

    /** parameters */
    parRelaxIntegrality_ = par_->getBoolPtrParam("RELAX_INTEGRALITY");
	gapTol_ = par_->getDblParam("DD/SUB/GAPTOL");

	/** augmented subproblem index */
    augs[0] = sind_;

    /** for auxiliary term */
    clbd_aux[0] = -COIN_DBL_MAX;
    cubd_aux[0] = COIN_DBL_MAX;
    obj_aux[0] = 0.0;

    /** decompose model */
    DSP_RTN_CHECK_THROW(
            model_->decompose(1, augs, 1, clbd_aux, cubd_aux, obj_aux,
                              mat, clbd, cubd, ctype, obj, rlbd, rubd));

    DSP_RTN_CHECK_THROW(
            model_->decomposeCoupling(1, augs, cpl_mat_, cpl_cols_, cpl_ncols));

	/** keep the original objective coefficient */
	obj_ = new double [mat->getNumCols()];
	CoinCopyN(obj, mat->getNumCols(), obj_);
	DSPdebugMessage("mat->getNumCols() = %d\n", mat->getNumCols());

    /** number of coupling variables and constraints for this subproblem */
    ncols_coupling_ = model_->getNumSubproblemCouplingCols(sind_);
    nrows_coupling_ = model_->getNumSubproblemCouplingRows(sind_);
    assert(ncols_coupling_ == cpl_ncols);
    assert(nrows_coupling_ == cpl_mat_->getNumRows());

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
			obj[j] *= probability;
		if (model_->isDro()) {
			if (sind_ < tssModel->getNumReferences()) {
				for (int j = tssModel->getNumCols(0); j < tssModel->getNumCols(0) + tssModel->getNumCols(1); ++j) {
					obj[j] *= tssModel->getReferenceProbability(sind_) / probability;
				}
			} else {
				CoinZeroN(obj + tssModel->getNumCols(0), tssModel->getNumCols(1));
			}
		}
		for (int j = 0; j < tssModel->getNumCols(1); ++j)
			obj_[tssModel->getNumCols(0)+j] /= probability;

#ifdef DSP_DEBUG
		DSPdebugMessage("sind_ = %d, probability = %e, lambdas = \n", sind_, model_->isDro() ? tssModel->getReferenceProbability(sind_) : probability);
		DspMessage::printArray(tssModel->getNumCols(0), obj);
#endif

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

	switch (par_->getIntParam("DD/SUB/SOLVER")) {
	case OsiCpx:
#ifdef DSP_HAS_CPX
		osi_ = new DspOsiCpx();
#else
		throw CoinError("Cplex is not available.", "createProblem", "DdSub");
#endif
		break;
	case OsiGrb:
#ifdef DSP_HAS_GRB
		osi_ = new DspOsiGrb();
#else
		throw CoinError("Gurobi is not available.", "createProblem", "DdSub");
#endif
		break;
	case OsiScip:
#ifdef DSP_HAS_SCIP
		osi_ = new DspOsiScip();
#else
		throw CoinError("Scip is not available.", "createProblem", "DdSub");
#endif
		break;
	default:
		throw CoinError("Invalid parameter value", "createProblem", "DdSub");
		break;
	}

	/** set number of cores */
	osi_->setNumCores(par_->getIntParam("DD/SUB/THREADS"));

	/** set display */
    osi_->setLogLevel(par_->getIntParam("DD/SUB/SOLVER/LOG_LEVEL"));

    /** load problem */
#ifdef DSP_DEBUG
	printf("create subproblem (s=%d):\n", sind_);
	DspMessage::printArray(mat->getNumCols(), obj);
#endif
    getSiPtr()->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	for (int j = 0; j < mat->getNumCols(); ++j) {
		if (ctype[j] != 'C')
			getSiPtr()->setInteger(j);
	}
    DSPdebug(mat->verifyMtx(4));

    /** set solution gap tolerance */
	if (nIntegers > 0)
	    osi_->setRelMipGap(gapTol_);

    END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

    FREE_MEMORY

    return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** add cut generator */
DSP_RTN_CODE DdSub::addCutGenerator() {
	BGN_TRY_CATCH
#ifdef DSP_HAS_SCIP
	if (par_->getIntParam("DD/SUB/SOLVER") == OsiScip) {
		DspOsiScip * osiscip = dynamic_cast<DspOsiScip*>(osi_);
		OsiScipSolverInterface *si = osiscip->scip_;
		if (si) {
			/** create constraint handler */
			SCIPconshdlrBendersDd *conshdlr = new SCIPconshdlrBendersDd(si->getScip());
			conshdlr->setOriginalVariables(si->getNumCols(), si->getScipVars(), 1);
			DSPdebugMessage("numcols %d, conshdlr %p\n", si->getNumCols(), conshdlr);
			/** add constraint handler */
			SCIP_CALL_ABORT(SCIPincludeObjConshdlr(si->getScip(), conshdlr, true));
			/* create constraint */
			SCIP_CONS * cons = NULL;
			SCIP_CALL_ABORT(SCIPcreateConsBenders(si->getScip(), &cons, "BendersDd"));
			SCIP_CALL_ABORT(SCIPaddCons(si->getScip(), cons));
			SCIP_CALL_ABORT(SCIPreleaseCons(si->getScip(), &cons));
		}
	}
#endif
	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)

    return DSP_RTN_OK;
}

/** update problem */
DSP_RTN_CODE DdSub::updateProblem(
		double * lambda,
		double probability,
		double primal_bound)
{
#define FREE_MEMORY       \
	FREE_ARRAY_PTR(newobj);

	double * newobj = NULL;

	BGN_TRY_CATCH

	int ncols = getSiPtr()->getNumCols();

	if (lambda)
	{
		/** allocate memory */
		newobj = new double [ncols];
		obj_offset_ = 0;

		/** update objective coefficients */
		assert(obj_);
		if (model_->isStochastic()) {
			// coefficients for the first-stage variables
			assert(nrows_coupling_<ncols);

			// coefficients for the second-stage variables
			CoinCopyN(lambda, nrows_coupling_, newobj);
			// printf("DdSub::updateProblem lambda:\n");
			// DspMessage::printArray(nrows_coupling_, lambda);
			for (int j = nrows_coupling_; j < ncols-1; ++j) {
				newobj[j] = obj_[j] * probability;
				// printf("update subproblem %d: obj_[%d] = %e, probability = %e, newobj -> %e\n", sind_, j, obj_[j], probability, newobj[j]);
			}
			newobj[ncols-1] = obj_[ncols-1];
		} else {
			CoinCopyN(obj_, ncols, newobj);

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
		}
#ifdef DSP_DEBUG
		printf("### newobj[%d]:\n", sind_);
		DspMessage::printArray(ncols, newobj);
#endif
		getSiPtr()->setObjective(newobj);
	}

	/** update primal bound (bounds of auxiliary constraint) */
	if (primal_bound < COIN_DBL_MAX)
		getSiPtr()->setColUpper(ncols - 1, primal_bound);//getSiPtr()->setColBounds(ncols - 1, primal_bound, primal_bound);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** push cuts */
DSP_RTN_CODE DdSub::pushCuts(OsiCuts * cuts)
{
	return DSP_RTN_OK;
}

/** set wall clock time limit */
void DdSub::setTimeLimit(double sec)
{
	assert(osi_);
	osi_->setTimeLimit(sec);
}

/** set accuracy tolerance */
void DdSub::setGapTol(double tol)
{
	assert(osi_);
	gapTol_ = tol;
	osi_->setRelMipGap(gapTol_);
}

/** set print level */
void DdSub::setPrintLevel(int level)
{
	osi_->setLogLevel(par_->getIntParam("DD/SUB/SOLVER/LOG_LEVEL"));
}
