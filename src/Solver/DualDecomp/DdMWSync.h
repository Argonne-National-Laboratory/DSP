/*
 * DdMWSync.h
 *
 *  Created on: Mar 23, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDMWSYNC_H_
#define SRC_SOLVER_DUALDECOMP_DDMWSYNC_H_

#include <Solver/DualDecomp/DdMWPara.h>

class DdMWSync: public DdMWPara {

public:

	/** constructor */
	DdMWSync(
			MPI_Comm     comm,   /**< MPI communicator */
			DecModel *   model,  /**< model pointer */
			DspParams *  par,    /**< parameters */
			DspMessage * message /**< message pointer */);

	/** destructor */
	virtual ~DdMWSync();

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** finalize */
	virtual DSP_RTN_CODE finalize();

protected:

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker();

private:

	/** sync Benders cut info */
	DSP_RTN_CODE syncBendersInfo(
			Solutions solutions, /**< solutions at which cuts are generated */
			OsiCuts & cuts       /**< cuts generated */);

	/** sync upper bound */
	DSP_RTN_CODE syncUpperbound(
			int nsolutions = 0, /**< number of solutions */
			vector<double> upperbounds = vector<double>() /**< list of upper bounds */);

	/** calculate upper bound */
	DSP_RTN_CODE calculateUpperbound(
			Solutions solutions, /**< solutions to evaluate */
			vector<double>&  upperbounds /**< list of upper bounds */);

	/** set coupling solutions */
	DSP_RTN_CODE setCouplingSolutions(
			Solutions &solutions /**< solution placeholder */);

	/** broadcast coupling solutions */
	DSP_RTN_CODE bcastCouplingSolutions(
			Solutions & solutions /**< solutions to broadcast */);

	/** scatter coupling solution */
	DSP_RTN_CODE scatterCouplingSolutions(
			Solutions & solutions /**< received solution placeholder */);
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWSYNC_H_ */
