/*
 * SCIPconshdlrBendersDd.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include <Utility/DspMessage.h>
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"

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
	clearCuts(cutsToAdd_);
	clearCuts(cutsAdded_);

	nvars_ = 0;

	return SCIP_OKAY;
}

/** set cutsToAdd_ */
void SCIPconshdlrBendersDd::setCutsToAdd(OsiCuts * cuts)
{
	/** cleaning */
	clearCuts(cutsToAdd_);
	for (int i = 0; i < cuts->sizeCuts(); ++i)
	{
		/** deep copy */
		cutsToAdd_->insert(cuts->rowCut(i));
	}
}

/** clean cutsAdded */
void SCIPconshdlrBendersDd::clearCuts(OsiCuts * cuts)
{
	for (int i = 0; i < cuts->sizeCuts(); ++i)
	{
		OsiRowCut * rc = cuts->rowCutPtr(i);
		FREE_PTR(rc);
	}
	cuts->dumpCuts();
}

void SCIPconshdlrBendersDd::generateCuts(
	int size,  /**< [in] size of x */
	double *x, /**< [in] master solution */
	OsiCuts *cuts /**< [out] cuts generated */)
{
	if (cutsToAdd_->sizeCuts() <= 0) return;

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
		if (row.getNumElements() == 0) continue;

		/** is optimality cut? */
		DSPdebug(rc->print());
		DSPdebugMessage("row.getNumElements() %d\n", row.getNumElements());
		bool isOptimalityCut = row.getIndices()[row.getNumElements() - 1] == nvars_ - 1;

		/** calculate efficacy */
//		DspMessage::printArray(size, x);
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
		cuts->insert(cutsToAdd_->rowCut(maxEfficaciousCut));

		/** mark as cut from pool */
		isCutFromPool = true;

		/** increment counter */
		nCutsReused_++;
	}

	END_TRY_CATCH(;)
}
