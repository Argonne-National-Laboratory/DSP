/*
 * DspModel.h
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_TREESEARCH_DSPMODEL_H_
#define SRC_TREESEARCH_DSPMODEL_H_

/** Coin */
#include "AlpsModel.h"
/** Dsp */
#include "Solver/DecSolver.h"

/**
 * This implements a model class for Coin-Alps library. This class provides a wrapper
 * for model class and data specific to a decomposition method.
 */
class DspModel: public AlpsModel {
public:

	/** default constructor */
	DspModel();

	/** constructor with solver */
	DspModel(DecSolver* solver);

	/** default destructor */
	virtual ~DspModel();

	/** solve model */
    DSP_RTN_CODE solve();

    /** Create the root node. Default: do nothing */
    virtual AlpsTreeNode * createRoot();

    /** Return true if all nodes on this process can be fathomed.*/
    virtual bool fathomAllNodes();

public:

    /** returns a pointer to the solver object */
    DecSolver* getSolver() {return solver_;}

private:

    DecSolver* solver_; /**< decomposition solver */
};

#endif /* SRC_TREESEARCH_DSPMODEL_H_ */
