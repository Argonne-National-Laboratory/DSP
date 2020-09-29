/*
 * SCIPconshdlrBendersDd.h
 *
 *  Created on: Dec 14, 2014
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_SCIPCONSHDLRBENDERSDD_H_
#define SRC_SOLVER_SCIPCONSHDLRBENDERSDD_H_

/** SCIP */
#include "objscip/objconshdlr.h"

/** COIN */
#include "OsiCuts.hpp"

/** DSP */
#include "Solver/Benders/SCIPconshdlrBenders.h"

class SCIPconshdlrBendersDd : public SCIPconshdlrBenders
{
public:

	/** default constructor */
	SCIPconshdlrBendersDd(SCIP * scip):
		SCIPconshdlrBenders(scip, "BendersDd", -2000000), nCutsReused_(0)
	{
		cutsToAdd_ = new OsiCuts;
		cutsAdded_ = new OsiCuts;
	}

	/** copy constructor */
	SCIPconshdlrBendersDd(SCIPconshdlrBendersDd * rhs):
		SCIPconshdlrBenders(rhs->scip_, "BendersDd", -2000000)
	{
		nCutsReused_ = rhs->getNumCutsReused();
		cutsToAdd_ = new OsiCuts(*rhs->getCutsToAdd());
		cutsAdded_ = new OsiCuts(*rhs->getCutsAdded());
	}

	/** default constructor */
	virtual ~SCIPconshdlrBendersDd()
	{
		if (cutsToAdd_) {
			clearCuts(cutsToAdd_);
			delete cutsToAdd_;
		}
		if (cutsAdded_) {
			clearCuts(cutsAdded_);
			delete cutsAdded_;
		}
	}

	/** destructor of constraint handler to free user data (called when SCIP is exiting) */
	virtual SCIP_DECL_CONSFREE(scip_free);

public:

	/** get cutsAdded_ */
	const OsiCuts * getCutsAdded() {return cutsAdded_;}

	/** get cutsToAdd_ */
	const OsiCuts * getCutsToAdd() {return cutsToAdd_;}

	/** get number of cuts reused */
	int getNumCutsReused() {return nCutsReused_;}

	/** set cutsToAdd_ */
	void setCutsToAdd(OsiCuts * cuts);

	/** clean cuts */
	void clearCuts(OsiCuts * cuts);

protected:

	/** generate Benders cuts */
	virtual void generateCuts(
		int size,  /**< [in] size of x */
		double *x, /**< [in] master solution */
		OsiCuts *cuts /**< [out] cuts generated */);

private:

	OsiCuts * cutsToAdd_;   /**< cut pool */
	OsiCuts * cutsAdded_;   /**< cut pool */
	int       nCutsReused_; /**< number of cuts reused */

};

#endif /* SRC_SOLVER_SCIPCONSHDLRBENDERSDD_H_ */
