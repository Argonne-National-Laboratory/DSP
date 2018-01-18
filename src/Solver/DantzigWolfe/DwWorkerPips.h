/*
 * DwWorkerPips.h
 *
 *  Created on: Dec 15, 2017
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWWORKERPIPS_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWWORKERPIPS_H_

#include "Solver/DantzigWolfe/DwWorkerMpi.h"
#include "Solver/DantzigWolfe/PipsInterface.h"

class DwWorkerPips: public DwWorkerMpi {
public:

	/** default constructor */
	DwWorkerPips(
			DecModel * model,
			DspParams * par,
			DspMessage * message,
			MPI_Comm comm);

	/** default destructor */
	virtual ~DwWorkerPips() {
		delete pips_;
		pips_ = NULL;
	}

	/** In this function, non-root processes receive signals
	 * from the root process and do proper processing. */
	virtual DSP_RTN_CODE receiver();

	enum {
		sig_initPips = 100,
		sig_solvePips,
		sig_clearMats	
	};

	/** initialize PIPS */
	virtual void initPips(int nscen = 0, int ncols = 0);

	/** Solve the master together using PIPS */
	virtual DSP_RTN_CODE solvePips(double weight = 0.0);

	/** clear rows added to PIPS matrices */
	virtual void clearMats();

	PipsInterface* pips() { return pips_; }

protected:
	PipsInterface* pips_;
	std::vector<int> rstart_; /**< to count number of rows added for each scenario */
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWWORKERPIPS_H_ */
