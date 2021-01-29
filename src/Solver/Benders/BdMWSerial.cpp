/*
 * BdMWSerial.cpp
 *
 *  Created on: Jul 10, 2016
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "Solver/Benders/BdMWSerial.h"
#include "Solver/Benders/SCIPconshdlrBenders.h"
#include "Solver/Benders/SCIPconshdlrDrBenders.h"
#include "Solver/Benders/SCIPconshdlrIntBenders.h"
#include "SolverInterface/DspOsiScip.h"

BdMWSerial::BdMWSerial(
		DecModel *   model,  /**< model pointer */
		DspParams *  par,    /**< parameters */
		DspMessage * message /**< message pointer */):
BdMW(model, par, message)
{
}

BdMWSerial::~BdMWSerial() {
	// TODO Auto-generated destructor stub
}

DSP_RTN_CODE BdMWSerial::init()
{
	BGN_TRY_CATCH

	/** set parameters */
	par_->setIntPtrParamSize("ARR_PROC_IDX", model_->getNumSubproblems());
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
		par_->setIntPtrParam("ARR_PROC_IDX", s, s);

	/** create master */
	master_ = new BdMaster(model_, par_, message_);

	/** initialize master */
	master_->init();

	/** create Benders worker */
	worker_ = new BdWorker(model_, par_, message_);

	if (master_->is_binary() == false && worker_->getBdSubPtr()->has_integer())
		warning_relaxation();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMWSerial::run()
{
	BGN_TRY_CATCH

	/** set constraint handler */
	
	if (par_->getIntParam("BD/MASTER/SOLVER")==OsiScip){
		DSPdebugMessage("setConshdlr for SCIP\n");
		master_->setConshdlr(constraintHandler());
	}
	else if (par_->getIntParam("BD/MASTER/SOLVER")==OsiGrb){
		master_->setBendersCallback(BendersCallbackFunc());
	}
	

	/** solve */
	DSPdebugMessage("solve\n");
	master_->solve();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE BdMWSerial::finalize()
{
	BGN_TRY_CATCH

	/** finalize and free master */
	master_->finalize();
	FREE_PTR(master_);

	/** free worker */
	FREE_PTR(worker_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

SCIPconshdlrBenders* BdMWSerial::constraintHandler()
{
	SCIPconshdlrBenders * conshdlr = NULL;

	BGN_TRY_CATCH

	int naux = par_->getIntParam("BD/NUM_CUTS_PER_ITER");
	int priority = par_->getIntParam("BD/CUT_PRIORITY");
	int sepa_solver = par_->getIntParam("BD/SUB/SOLVER");

	/** get solver interface */
	OsiScipSolverInterface * si = dynamic_cast<OsiScipSolverInterface*>(master_->getSiPtr());

	/** get Benders subproblem solver */
	BdSub *bdsub = worker_->getBdSubPtr();

	/** check whether integer Benders cuts need used or not */
	bool add_integer_benders = false;
	if (model_->isStochastic() && master_->is_binary())
		add_integer_benders = true;

	/** Benders constraint handler */
	if (add_integer_benders){
		// Distributionally robust integer benders is implemented in SCIPconshdlrIntBenders.
		conshdlr = new SCIPconshdlrIntBenders(si->getScip(), "Benders", priority, sepa_solver);
	}
	else
	{
		if (model_->isDro())
			conshdlr = new SCIPconshdlrDrBenders(si->getScip(), "Benders", priority, sepa_solver);
		else{
			si->getScip();
			conshdlr = new SCIPconshdlrBenders(si->getScip(), "Benders", priority);
		}
	}
	conshdlr->setDecModel(model_);
	conshdlr->setBdSub(bdsub);
	conshdlr->setOriginalVariables(si->getNumCols(), si->getScipVars(), naux);

	END_TRY_CATCH_RTN(;, NULL)

	return conshdlr;
}

DSP_RTN_CODE BdMWSerial::BendersCallbackFunc(){

	BGN_TRY_CATCH

	BendersCallback * Bdcb = NULL;

	BGN_TRY_CATCH

	int naux = par_->getIntParam("BD/NUM_CUTS_PER_ITER");
	int priority = par_->getIntParam("BD/CUT_PRIORITY");
	int sepa_solver = par_->getIntParam("BD/SUB/SOLVER");

	/** get solver interface */
	DspOsi * osi = dynamic_cast<DspOsi*>(master_->getDspOsiPtr());

	/** get Benders subproblem solver */
	BdSub *bdsub = worker_->getBdSubPtr();

	/** check whether integer Benders cuts need used or not */
	bool add_integer_benders = false;
	if (model_->isStochastic() && master_->is_binary())
		add_integer_benders = true;

	/** Benders constraint handler 
	 * TODO: add integer Benders cuts
	 * TODO: if Dro
	*/
	
	if (model_->isDro())
		printf("To implement Benders Callback for Dro...\n");
	else{
		Bdcb = new BendersCallback(osi, "Benders", priority);
	}
	
	Bdcb->setDecModel(model_);
	Bdcb->setBdSub(bdsub);
	Bdcb->setOriginalVariables(osi->si_->getNumCols(), naux);

	END_TRY_CATCH_RTN(;, NULL)

	return conshdlr;

	return DSP_RTN_OK
}