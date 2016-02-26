/*
 * DdSolver.h
 *
 *  Created on: Feb 9, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DUALDECOMP_DDSOLVER_H_
#define SRC_SOLVER_DUALDECOMP_DDSOLVER_H_

#include "Solver/DecSolver.h"

class DdSolver: public DecSolver {
public:

	/** constructor */
	DdSolver(DspParams * par, DecModel * model, StoMessage * message):
		DecSolver(par, model, message),
		scount_(0), rcount_(0), sbuf_(NULL), rbuf_(NULL) {}

	/** destructor */
	virtual ~DdSolver()
	{
		FREE_ARRAY_PTR(sbuf_);
		FREE_ARRAY_PTR(rbuf_);
	}

	/** solve */
	virtual STO_RTN_CODE solve() {return STO_RTN_OK;}

	/** receive message from worker */
	virtual STO_RTN_CODE recvMessage(int source, int size, double * message) {return STO_RTN_OK;}

	/** get size of receiving message from worker */
	virtual int getRecvCount() {return rcount_;}

	/** get size of sending message from worker */
	virtual int getSendCount() {return scount_;}

	/** get sending message from worker */
	virtual const double * getSendMessage() {return sbuf_;}

	/** initialize sending message */
	virtual STO_RTN_CODE initSendMessage(int size, double * buf = NULL)
	{
		BGN_TRY_CATCH

		FREE_ARRAY_PTR(sbuf_);
		scount_ = size;
		sbuf_ = new double [scount_];
		if (buf != NULL)
			CoinCopyN(buf, scount_, sbuf_);

		END_TRY_CATCH_RTN(;,STO_RTN_ERR)

		return STO_RTN_OK;
	}

	/** initialize receiving message */
	virtual STO_RTN_CODE initRecvMessage(int size, double * buf = NULL)
	{
		BGN_TRY_CATCH

		FREE_ARRAY_PTR(rbuf_);
		rcount_ = size;
		rbuf_ = new double [rcount_];
		if (buf != NULL)
			CoinCopyN(buf, rcount_, rbuf_);

		END_TRY_CATCH_RTN(;,STO_RTN_ERR)

		return STO_RTN_OK;
	}

protected:

	/** create problem */
	virtual STO_RTN_CODE createProblem() {return STO_RTN_OK;}

	/** create message */
	virtual STO_RTN_CODE createMessage() {return STO_RTN_OK;}

protected:

	int scount_;    /**< send count */
	int rcount_;    /**< recv count */

	double * sbuf_; /**< send buffer */
	double * rbuf_; /**< recv buffer */
};

#endif /* SRC_SOLVER_DUALDECOMP_DDSOLVER_H_ */
