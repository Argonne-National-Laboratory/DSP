/*
 * SCIPconshdlrBendersDd.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#include "Utility/StoMessage.h"

/** constraint data for Benders cuts */
struct SCIP_ConsData
{
	/** nothing to do */
};

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
SCIP_DECL_CONSFREE(SCIPconshdlrBendersDd::scip_free)
{
	DSPdebugMessage("scip_free\n");

	/** release variables */
	for (int j = 0; j < nvars_; ++j)
		vars_[j] = NULL;
	SCIPfreeMemoryArray(scip_, &vars_);

	/** release cut pool */
	cutsToAdd_->dumpCuts();
	for (int i = 0; i < cutsAdded_->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cutsAdded_->rowCutPtr(i);
		FREE_PTR(rc);
	}
	cutsAdded_->dumpCuts();

	nvars_ = 0;

	return SCIP_OKAY;
}

/** clone method which will be used to copy constraint handler and variable pricer objects */
SCIP_DECL_CONSHDLRCLONE(scip::ObjProbCloneable* SCIPconshdlrBendersDd::clone)
{
	*valid = true;
	SCIPconshdlrBendersDd * conshdlrclone = new SCIPconshdlrBendersDd(const_cast<SCIPconshdlrBendersDd*>(this));
	conshdlrclone->setOriginalVariables(nvars_, vars_, 1);
	return conshdlrclone;
}

/** set cutsToAdd_ */
void SCIPconshdlrBendersDd::setCutsToAdd(OsiCuts * cuts)
{
	cutsToAdd_->dumpCuts();
	for (int i = 0; i < cuts->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts->rowCutPtr(i);
		if (!rc) continue;

		/** shallow copy */
		cutsToAdd_->insert(rc);
	}
}

/** clean cutsAdded */
void SCIPconshdlrBendersDd::clearCutsAdded()
{
	for (int i = 0; i < cutsAdded_->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cutsAdded_->rowCutPtr(i);
		FREE_PTR(rc);
	}
	cutsAdded_->dumpCuts();
}

void SCIPconshdlrBendersDd::generateCuts(
		int size,      /**< [in] size of x */
		double * x,    /**< [in] master solution */
		int where,     /**< [in] where to be called */
		OsiCuts * cuts /**< [out] cuts generated */)
{
	BGN_TRY_CATCH

	/** update cut efficacy in cutsToAdd_ */
	bool isCutFromPool = false;
	int maxEfficaciousCut = -1;
	double maxEfficacy = 1.e-6;
	DSPdebugMessage("Number of cutsToAdd_ %d\n", cutsToAdd_->sizeCuts());
	for (int i = cutsToAdd_->sizeCuts() - 1; i >= 0; --i)
	{
		OsiRowCut * rc = cutsToAdd_->rowCutPtr(i);
		if (!rc) continue;
		const CoinPackedVector row = rc->row();

		/** is optimality cut? */
		DSPdebugMessage("row.getNumElements() %d\n", row.getNumElements());
		bool isOptimalityCut = row.getIndices()[row.getNumElements() - 1] == nvars_ - 1;

		/** calculate efficacy */
		double efficacy = rc->violated(x);
		if (isOptimalityCut) efficacy /= row.twoNorm();

		/** determine if efficacious */
		if (efficacy <= 1.e-6)
			continue;
		else if (efficacy > maxEfficacy)
		{
			maxEfficaciousCut = i;
			maxEfficacy = efficacy;
		}
	}

	DSPdebugMessage("maxEfficacy %e\n", maxEfficacy);
	if (maxEfficaciousCut >= 0)
	{
		/** move one from cut pool */
		OsiRowCut * rc = NULL;
		if (where == from_scip_check)
			rc = cutsToAdd_->rowCutPtr(maxEfficaciousCut);
		else
			rc = cutsToAdd_->rowCutPtrAndZap(maxEfficaciousCut);
		cuts->insert(rc);

		/** mark as cut from pool */
		isCutFromPool = true;

		/** increment counter */
		nCutsReused_++;
	}

	END_TRY_CATCH(;)
}
