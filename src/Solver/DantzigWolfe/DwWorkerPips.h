/*
 * DwWorkerPips.h
 *
 *  Created on: Dec 15, 2017
 *      Author: Kibaek Kim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWWORKERPIPS_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWWORKERPIPS_H_

#include "Model/TssModel.h"
#include "Solver/DantzigWolfe/DwCol.h"
#include "Solver/DantzigWolfe/DwWorkerMpi.h"
#include "Solver/DantzigWolfe/DwBundlePipsInput.h"

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
		tss_ = NULL;
		delete pips_;
		pips_ = NULL;
	}

	/** In this function, non-root processes receive signals
	 * from the root process and do proper processing. */
	virtual DSP_RTN_CODE receiver();

	enum {
		sig_sync = 100,
		sig_solvePips	
	};

	/** run PIPS */
	virtual DSP_RTN_CODE sync(std::vector<double> bestdualsol, std::vector<DwCol*> dwcols);

	/** Solve the master together using PIPS */
	virtual DSP_RTN_CODE solvePips(double weight);

protected:
	TssModel* tss_;
	DwBundlePipsInput* pips_;

public:

	/** PIPS solutions */
	double pips_objval_;
	std::vector<std::vector<double>> pips_beta_;
	std::vector<std::vector<double>> pips_theta_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWWORKERPIPS_H_ */
