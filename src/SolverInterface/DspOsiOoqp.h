/**
 * DspOsiOoqp.h
 *
 * 5/28/2020
 * Kibaek Kim
 */

#ifndef SRC_SOLVERINTERFACE_DSPOSIOOQP_H_
#define SRC_SOLVERINTERFACE_DSPOSIOOQP_H_

#ifdef DSP_HAS_OOQP
#include "SolverInterface/DspOsi.h"
#include "SolverInterface/OsiOoqpSolverInterface.hpp"

class DspOsiOoqp : public DspOsi {
public:

	/** default constructor */
	DspOsiOoqp() {
        si_ = new OsiOoqpSolverInterface();
        ooqp_ = dynamic_cast<OsiOoqpSolverInterface*>(si_);
    }

	/** copy constructor */
	DspOsiOoqp(const DspOsiOoqp& rhs) {
        si_ = new OsiOoqpSolverInterface(*(rhs.ooqp_));
        ooqp_ = dynamic_cast<OsiOoqpSolverInterface*>(si_);}

	/** clone constructor */
	virtual DspOsiOoqp* clone() const {
        return new DspOsiOoqp(*this);
    }

	/** destructor */
	virtual ~DspOsiOoqp() {
		ooqp_ = NULL;
	}

	/** load quadratic objective */
	virtual void loadQuadraticObjective(const CoinPackedMatrix &mat) {
		throw CoinError("Quadratic objective is not supported.", "loadQuadraticObjective", "DspOsi");
	}

	/** solve problem */
	virtual void solve() {
        si_->resolve();
    }

	virtual void use_barrier() {}

	OsiOoqpSolverInterface *ooqp_;
};

#endif

#endif /* SRC_SOLVERINTERFACE_DSPOSIOOQP_H_ */