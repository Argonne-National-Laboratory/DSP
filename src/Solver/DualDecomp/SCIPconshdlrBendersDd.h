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
		SCIPconshdlrBenders(scip, -2000000), nCutsReused_(0)
	{
		cutsToAdd_ = new OsiCuts;
		cutsAdded_ = new OsiCuts;
	}

	/** copy constructor */
	SCIPconshdlrBendersDd(SCIPconshdlrBendersDd * rhs):
		SCIPconshdlrBenders(rhs->scip_, -2000000)
	{
		nCutsReused_ = rhs->getNumCutsReused();
		cutsToAdd_ = new OsiCuts(*rhs->getCutsToAdd());
		cutsAdded_ = new OsiCuts(*rhs->getCutsAdded());
	}

	/** default constructor */
	virtual ~SCIPconshdlrBendersDd()
	{
		/** nothing to do */
	}

	/** destructor of constraint handler to free user data (called when SCIP is exiting) */
	virtual SCIP_DECL_CONSFREE(scip_free);

	/** clone method which will be used to copy constraint handler and variable pricer objects */
	virtual SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable* clone);

public:

	/** get cutsAdded_ */
	const OsiCuts * getCutsAdded() {return cutsAdded_;}

	/** get cutsToAdd_ */
	const OsiCuts * getCutsToAdd() {return cutsToAdd_;}

	/** get number of cuts reused */
	int getNumCutsReused() {return nCutsReused_;}

	/** set cutsToAdd_ */
	void setCutsToAdd(OsiCuts * cuts);

	/** clean cutsAdded */
	void clearCutsAdded();

protected:

	/** generate Benders cuts */
	virtual void generateCuts(
			int size,      /**< [in] size of x */
			double * x,    /**< [in] master solution */
			int where,     /**< [in] where to be called */
			OsiCuts * cuts /**< [out] cuts generated */);

private:

	OsiCuts * cutsToAdd_;   /**< cut pool */
	OsiCuts * cutsAdded_;   /**< cut pool */
	int       nCutsReused_; /**< number of cuts reused */

};

#endif /* SRC_SOLVER_SCIPCONSHDLRBENDERSDD_H_ */
