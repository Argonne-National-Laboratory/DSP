/*
 * DdSub.h
 *
 *  Renamed on: Feb 16, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/DdSub.h"
#include "Utility/StoMacros.h"
#include "Utility/StoMessage.h"
#include "SolverInterface/SolverInterfaceClp.h"
#include "SolverInterface/SolverInterfaceScip.h"
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#include "SolverInterface/SCIPbranchruleLB.h"

DdSub::~DdSub()
{
	FREE_PTR(si_);
	FREE_PTR(obj_);
	FREE_ARRAY_PTR(lambda_);
	FREE_ARRAY_PTR(cpl_rhs_);
	FREE_2D_ARRAY_PTR(nsols_, solutions_);
//	for (unsigned i = 0; i < s_primsols_.size(); ++i)
//		FREE_ARRAY_PTR(s_primsols_[i]);
	nsols_ = 0;
	obj_offset_ = 0;
	parRelaxIntegrality_ = NULL;
}

/** create problem */
STO_RTN_CODE DdSub::createProblem(DecModel * model)
{
#define FREE_MEMORY       \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	/** problem data */
	CoinPackedMatrix * mat = NULL;
	double * clbd = NULL;
	double * cubd = NULL;
	char * ctype  = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;
	int augs[1];
	double clbd_aux[1];
	double cubd_aux[1];
	double obj_aux[1];
	int cpl_ncols;

	BGN_TRY_CATCH

	/** parameters */
	parRelaxIntegrality_ = par_->getBoolPtrParam("RELAX_INTEGRALITY");
	gapTol_              = par_->getDblParam("SCIP/GAP_TOL");

	/** augmented subproblem index */
	augs[0] = sind_;

	/** for auxiliary term */
	clbd_aux[0] = -COIN_DBL_MAX;
	cubd_aux[0] = COIN_DBL_MAX;
	obj_aux[0] = 0.0;

	/** decompose model */
	STO_RTN_CHECK_THROW(
			model->decompose(1, augs, 1, clbd_aux, cubd_aux, obj_aux,
					mat, clbd, cubd, ctype, obj_, rlbd, rubd),
			"decompose", "DecModel");

	STO_RTN_CHECK_THROW(
			model->decomposeCoupling(1, augs, cpl_mat_, cpl_cols_, cpl_ncols),
			"decomposeCoupling", "DecModel");

	/** number of coupling variables and constraints for this subproblem */
	ncols_coupling_ = model->getNumSubproblemCouplingCols(sind_);
	nrows_coupling_ = model->getNumSubproblemCouplingRows(sind_);
	assert(ncols_coupling_ == cpl_ncols);
	assert(nrows_coupling_ == cpl_mat_->getNumRows());

	/** storage for lambda */
	lambda_ = new double [nrows_coupling_];
	CoinZeroN(lambda_, nrows_coupling_);

	/** copy right-hand side of coupling rows */
	cpl_rhs_ = new double [nrows_coupling_];
	for (int i = 0; i < nrows_coupling_; i++)
		cpl_rhs_[i] = model->getRhsCouplingRow(i);
	obj_offset_ = 0;

	/** number of integer variables */
	int nIntegers = model->getNumIntegers();

	/** adjust first-stage cost */
	if (model->isStochastic())
	{
		TssModel * tssModel;
		try
		{
			tssModel = dynamic_cast<TssModel *>(model);
		}
		catch (const std::bad_cast& e)
		{
			printf("Error: Model claims to be stochastic when it is not\n");
			return STO_RTN_ERR;
		}

		double probability = tssModel->getProbability()[sind_];
		for (int j = 0; j < tssModel->getNumCols(0); ++j)
			obj_[j] *= probability;

		/** convert column types */
		if (parRelaxIntegrality_[0])
		{
			for (int j = 0; j < tssModel->getNumCols(0); ++j)
			{
				if (ctype[j] != 'C')
					nIntegers--;
				ctype[j] = 'C';
			}
		}
		if (parRelaxIntegrality_[1])
		{
			for (int j = 0; j < tssModel->getNumCols(1); ++j)
			{
				if (ctype[tssModel->getNumCols(0) + j] != 'C')
					nIntegers--;
			}
			CoinFillN(ctype + tssModel->getNumCols(0), tssModel->getNumCols(1), 'C');
		}
	}
	else
	{
		if (parRelaxIntegrality_[0] || parRelaxIntegrality_[1])
		{
			for (int j = 0; j < mat->getNumCols(); j++)
			{
				if (ctype[j] != 'C')
					nIntegers--;
				ctype[j] = 'C';
			}
		}
	}

	if (nIntegers > 0)
		si_ = new SolverInterfaceScip(par_);
	else
		si_ = new SolverInterfaceClp(par_);

	/** load problem */
	si_->loadProblem(mat, clbd, cubd, obj_, ctype, rlbd, rubd, "DdSub");
	DSPdebug(mat->verifyMtx(4));

	/** set solution gap tolerance */
	si_->setGapTol(gapTol_);

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** add cut generator */
STO_RTN_CODE DdSub::addCutGenerator(BdSub * bdsub)
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		SCIPconshdlrBendersDd * conshdlr = new SCIPconshdlrBendersDd(si->getSCIP());
		conshdlr->setOriginalVariables(si->getNumCols(), si->getSCIPvars(), 1);

		/** add constraint handler */
		si->addConstraintHandler(conshdlr, true, true);
	}
	else
	{
		/** TODO */
		printf("Warning: Cut generation supports only SCIP.\n");
	}
	return STO_RTN_OK;
}

/** add branch rule */
STO_RTN_CODE DdSub::addBranchrule()
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		/** add branch rule */
		si->addBranchrule(new SCIPbranchruleLB(si->getSCIP()), true);
	}
	else
	{
		/** TODO */
		printf("Warning: Branch rule supports only SCIP.\n");
	}
	return STO_RTN_OK;
}

/** change branch rule */
STO_RTN_CODE DdSub::chgBranchrule(double lb)
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		SCIPbranchruleLB * branchrule = dynamic_cast<SCIPbranchruleLB*>(si->findObjBranchrule("knownLB"));
		if (branchrule)
		{
			branchrule->setLowerBound(lb);
		}
	}
	return STO_RTN_OK;
}

/** update problem */
STO_RTN_CODE DdSub::updateProblem(
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

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** solve problem */
STO_RTN_CODE DdSub::solve()
{
	si_->solve();

	return STO_RTN_OK;
}

/** free solution process data */
STO_RTN_CODE DdSub::freeSolve(bool restart)
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si) si->freeSolve(restart);

	return STO_RTN_OK;
}

/** free all solution process data */
STO_RTN_CODE DdSub::freeTransform()
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si) si->freeTransform();
	return STO_RTN_OK;
}

/** collect cuts */
STO_RTN_CODE DdSub::collectCuts(OsiCuts * cuts)
{
	assert(cuts);
	/** TODO only for SCIP */
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		const OsiCuts * cutsAdded = si->getCuts();
		assert(cutsAdded);

		DSPdebugMessage("cutsAdded %d\n", cutsAdded->sizeCuts());
		for (int i = 0; i < cutsAdded->sizeCuts(); ++i)
		{
			OsiRowCut rc = cutsAdded->rowCut(i);
			//rc.print();
			cuts->insertIfNotDuplicate(rc);
		}

		/** clear cuts added */
		si->clearCuts();
	}
	return STO_RTN_OK;
}

/** push cuts */
STO_RTN_CODE DdSub::pushCuts(OsiCuts * cuts)
{
	/** TODO only for SCIP */
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		si->setCuts(cuts);
	}
	return STO_RTN_OK;
}

/** collect solutions */
STO_RTN_CODE DdSub::collectSolutions()
{
	/** free memeory */
	FREE_2D_ARRAY_PTR(nsols_, solutions_);

	/** TODO only for SCIP */
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		/** get solutions */
		nsols_ = si->getSolutions(&solutions_);
	}
	else
	{
		nsols_ = 1;
		solutions_ = new double * [nsols_];
		CoinCopyN(si->getSolution(), si->getNumCols(), solutions_[0]);
	}

	return STO_RTN_OK;
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
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
		si_->setPrintLevel(level);
	else
		si_->setPrintLevel(CoinMax(0, level - 2));
}

/** get MPI message buffer: [scenario index, objval, coupling solution] */
double * DdSub::MPImsgbuf()
{
	double * msgbuf = new double [ncols_coupling_ + 2];
	MPImsgbuf(msgbuf);
	return msgbuf;
}

/** get MPI message buffer */
STO_RTN_CODE DdSub::MPImsgbuf(double * msgbuf)
{
	/** The first element in message buffer should be scenario index. */
	msgbuf[0] = static_cast<double>(sind_);

	/** The second element should be primal objective value. */
	msgbuf[1] = si_->getPrimalBound() + obj_offset_;

	/** The second element should be dual objective value. */
	msgbuf[2] = si_->getDualBound() + obj_offset_;

	/** The following elements are for the coupling solution. */
	const double * solution = si_->getSolution();
	if (solution != NULL)
	{
		for (int i = 0; i < ncols_coupling_; i++)
			msgbuf[3+i] = solution[cpl_cols_[i]];
		solution = NULL;
	}

	return STO_RTN_OK;
}
