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

	typedef vector<CoinPackedVector*> Solutions;

public:

	/** constructor */
	DdMWSync(
			MPI_Comm          comm,   /**< MPI communicator */
			DdMaster *        master, /**< master problem */
			vector<DdWorker*> worker  /**< worker for finding lower bounds */);

	/** destructor */
	virtual ~DdMWSync();

protected:

	/** initialize */
	virtual DSP_RTN_CODE init();

	/** run master process */
	virtual DSP_RTN_CODE runMaster();

	/** run worker processes */
	virtual DSP_RTN_CODE runWorker();

private:

	/** sync Benders cut info */
	DSP_RTN_CODE syncBendersInfo(
			Solutions solutions, /**< solutions at which cuts are generated */
			OsiCuts & cuts       /**< cuts generated */);

	/** generate Benders cuts */
	DSP_RTN_CODE generateBendersCuts(
			Solutions solutions, /**< solutions at which cuts are generated */
			OsiCuts & cuts       /**< cuts generated */);

	/** receive Benders cuts */
	DSP_RTN_CODE recvBendersCuts(OsiCuts & cuts);

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

	/** store coupling solution */
	DSP_RTN_CODE storeCouplingSolutions(Solutions & stored);

	/** send coupling solution */
	DSP_RTN_CODE sendCouplingSolutions(Solutions solutions);

	/** get coupling solutions */
	DSP_RTN_CODE getCouplingSolutions(
			Solutions &solutions /**< solution placeholder */);

	/** receive coupling solutions */
	DSP_RTN_CODE recvCouplingSolutions(
			Solutions &solutions /**< received solution placeholder */);

	/** scatter coupling solution */
	DSP_RTN_CODE scatterCouplingSolutions(
			Solutions & solutions /**< received solution placeholder */);
};

#endif /* SRC_SOLVER_DUALDECOMP_DDMWSYNC_H_ */
