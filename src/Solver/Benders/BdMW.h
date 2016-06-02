/*
 * BdMW.h
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_BENDERS_BDMW_H_
#define SRC_SOLVER_BENDERS_BDMW_H_

#include "Solver/BaseMasterWorker.h"
#include "Solver/Benders/BdMaster.h"
#include "Solver/Benders/BdWorker.h"

class BdMW: public BaseMasterWorker {
public:

	enum
	{
		MASTER_NEEDS_CUTS = 0,
		MASTER_STOPPED
	};

	/** constructor */
	BdMW(MPI_Comm comm, BdMaster * master, BdWorker * worker);

	/** destructor */
	virtual ~BdMW();

	/** run the framework */
	virtual DSP_RTN_CODE run();

protected:

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

public:

	/** get priaml solution */
	const double * getPrimalSolution() {return primsol_;}

protected:

	/** Common member variables */
	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	BdMaster * master_; /**< master */
	BdWorker * worker_; /**< worker */

private:

	double * primsol_;
};

#endif /* SRC_SOLVER_BENDERS_BDMW_H_ */
