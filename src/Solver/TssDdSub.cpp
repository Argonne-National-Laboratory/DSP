/*
 * TssDdSub.cpp
 *
 *  Created on: Dec 10, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/TssDdSub.h"
#include "Utility/StoMacros.h"
#include "Utility/StoMessage.h"
#include "Solver/SolverInterfaceClp.h"
#include "Solver/SolverInterfaceScip.h"
#include "Solver/SCIPconshdlrBendersDd.h"
#include "Solver/SCIPbranchruleLB.h"

TssDdSub::~TssDdSub()
{
	FREE_PTR(si_);
	FREE_PTR(obj_);
	FREE_ARRAY_PTR(lambda_);
	FREE_2D_ARRAY_PTR(nsols_, solutions_);
	nsols_ = 0;
}

/** create problem */
STO_RTN_CODE TssDdSub::createProblem(
		double     probability, /**< probability */
		TssModel * model)
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

	BGN_TRY_CATCH

	/** augmented scenario index */
	augs[0] = sind_;

	/** for auxiliary term */
	clbd_aux[0] = -COIN_DBL_MAX;
	cubd_aux[0] = COIN_DBL_MAX;
	obj_aux[0] = 0.0;

	/** decompose model */
	STO_RTN_CHECK_THROW(
			model->decompose(1, augs, 1, clbd_aux, cubd_aux, obj_aux,
					mat, clbd, cubd, ctype, obj_, rlbd, rubd), "decompose", "TssModel");

	/** number of first-stage variables */
	ncols_first_ = model->getNumCols(0);

	/** storage for lambda */
	lambda_ = new double [ncols_first_];
	CoinZeroN(lambda_, ncols_first_);

	/** adjust first-stage cost */
	for (int j = 0; j < ncols_first_; ++j)
		obj_[j] *= probability;

	/** number of integer variables in the core */
	int nIntegers = model->getNumCoreIntegers();

	/** convert column types */
	if (par_->relaxIntegrality_[0])
	{
		for (int j = 0; j < model->getNumCols(0); ++j)
		{
			if (ctype[j] != 'C')
				nIntegers--;
			ctype[j] = 'C';
		}
	}
	if (par_->relaxIntegrality_[1])
	{
		for (int j = 0; j < model->getNumCols(1); ++j)
		{
			if (ctype[model->getNumCols(0) + j] != 'C')
				nIntegers--;
		}
		CoinFillN(ctype + model->getNumCols(0), model->getNumCols(1), 'C');
	}

	if (nIntegers > 0)
		si_ = new SolverInterfaceScip(par_);
	else
		si_ = new SolverInterfaceClp(par_);

	/** load problem */
	si_->loadProblem(mat, clbd, cubd, obj_, ctype, rlbd, rubd, "TssDdSub");
	DSPdebug(mat->verifyMtx(4));

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** add cut generator */
STO_RTN_CODE TssDdSub::addCutGenerator(TssBdSub * tss)
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		SCIPconshdlrBendersDd * conshdlr = new SCIPconshdlrBendersDd(si->getSCIP());
		DSPdebugMessage("add cut generator %p\n", tss);
		conshdlr->assignData(tss, si->getNumCols(), ncols_first_, obj_);
		conshdlr->setOriginalVariables(si->getNumCols(), si->getSCIPvars());

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

/** change cut generator */
STO_RTN_CODE TssDdSub::chgCutGenerator(TssBdSub * tss)
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		SCIPconshdlrBendersDd * conshdlr = dynamic_cast<SCIPconshdlrBendersDd*>(si->findObjConshdlr("BendersDd"));
		if (conshdlr)
		{
			DSPdebugMessage("change cut generator to %p\n", tss);
			conshdlr->assignData(tss, si->getNumCols(), ncols_first_, obj_);
		}
	}
	return STO_RTN_OK;
}

/** add branch rule */
STO_RTN_CODE TssDdSub::addBranchrule()
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
STO_RTN_CODE TssDdSub::chgBranchrule(double lb)
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
STO_RTN_CODE TssDdSub::updateProblem(
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
		CoinCopyN(lambda, ncols_first_, lambda_);

		/** allocate memory */
		newobj = new double [ncols];

		/** update objective coefficients */
		assert(obj_);
		for (int j = 0; j < ncols_first_; ++j)
			newobj[j] = obj_[j] + lambda[j];
		for (int j = ncols_first_; j < ncols; ++j)
			newobj[j] = obj_[j];
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
STO_RTN_CODE TssDdSub::solve()
{
#if 0
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
	{
		si->enableUpperBoundingCuts(enableCuts, hasUb);
	}
#endif
	//si_->writeMps("subprob");
	si_->solve();

	return STO_RTN_OK;
}

/** free solution process data */
STO_RTN_CODE TssDdSub::freeSolve(bool restart)
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si) si->freeSolve(restart);

	return STO_RTN_OK;
}

/** free all solution process data */
STO_RTN_CODE TssDdSub::freeTransform()
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si) si->freeTransform();
	return STO_RTN_OK;
}

/** collect cuts */
STO_RTN_CODE TssDdSub::collectCuts(OsiCuts * cuts)
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
STO_RTN_CODE TssDdSub::pushCuts(OsiCuts * cuts)
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
STO_RTN_CODE TssDdSub::collectSolutions()
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
void TssDdSub::setTimeLimit(double sec)
{
	assert(si_);
	si_->setTimeLimit(sec);
}

/** set print level */
void TssDdSub::setPrintLevel(int level)
{
	SolverInterfaceScip * si = dynamic_cast<SolverInterfaceScip*>(si_);
	if (si)
		si_->setPrintLevel(level);
	else
		si_->setPrintLevel(CoinMax(0, level - 2));
}

/** get MPI message buffer: [scenario index, objval, first-stage solution] */
double * TssDdSub::MPImsgbuf()
{
	double * msgbuf = new double [ncols_first_ + 2];
	MPImsgbuf(msgbuf);
	return msgbuf;
}

/** get MPI message buffer */
STO_RTN_CODE TssDdSub::MPImsgbuf(double * msgbuf)
{
	/** The first element in message buffer should be scenario index. */
	msgbuf[0] = static_cast<double>(sind_);

	/** The second element should be objective value. */
	msgbuf[1] = si_->getPrimalBound();

	/** The following elements are for the first-stage solution. */
	const double * solution = si_->getSolution();
	if (solution != NULL)
	{
		CoinCopyN(solution, ncols_first_, msgbuf + 2);
		solution = NULL;
	}

	return STO_RTN_OK;
}

