/**
 * DspOsiOoqpEps.h
 *
 * 5/28/2020
 * Kibaek Kim
 */

#ifndef SRC_SOLVERINTERFACE_DSPOSIOOQPEPS_H_
#define SRC_SOLVERINTERFACE_DSPOSIOOQPEPS_H_

#ifdef DSP_HAS_OOQP
#include "SolverInterface/DspOsi.h"
#include "SolverInterface/OoqpEps.h"

class DspOsiOoqpEps : public DspOsi {
public:

	/** default constructor */
	DspOsiOoqpEps() {
        si_ = new OoqpEps();
        ooqp_ = dynamic_cast<OoqpEps*>(si_);
    }

	/** copy constructor */
	DspOsiOoqpEps(const DspOsiOoqpEps& rhs) {
        si_ = new OoqpEps(*(rhs.ooqp_));
        ooqp_ = dynamic_cast<OoqpEps*>(si_);}

	/** clone constructor */
	virtual DspOsiOoqpEps* clone() const {
        return new DspOsiOoqpEps(*this);
    }

	/** destructor */
	virtual ~DspOsiOoqpEps() {
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

	OoqpEps *ooqp_;
};

#endif

#endif /* SRC_SOLVERINTERFACE_DSPOSIOOQPEPS_H_ */