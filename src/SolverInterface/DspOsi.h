/**
 * DspOsi.h
 *
 * 12/12/2019
 * Kibaek Kim
 */

#ifndef SRC_SOLVERINTERFACE_DSPOSI_H_
#define SRC_SOLVERINTERFACE_DSPOSI_H_

#include "Utility/DspRtnCodes.h"
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGlpk.h"
#include "SolverInterface/DspOsiScip.h"

#ifdef DSP_HAS_OOQP
#include "SolverInterface/OsiOoqpSolverInterface.hpp"
#include "SolverInterface/OoqpEps.h"
#endif

/** convert Osi-status to dsp-status */
inline void convertOsiToDspStatus(OsiSolverInterface* si, int& status) {

#ifdef DSP_HAS_CPX
	OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si);
	if (cpx) {
		convertOsiCpxToDspStatus(cpx, status);
		return;
	}
#endif /* DSP_HAS_CPX */

#ifdef DSP_HAS_SCIP
	OsiScipSolverInterface* scip = dynamic_cast<OsiScipSolverInterface*>(si);
	if (scip) {
		convertOsiScipToDspStatus(scip, status);
		return;
	}
#endif /* DSP_HAS_SCIP */

	if (si->isProvenOptimal())
		status = DSP_STAT_OPTIMAL;
	else if (si->isProvenPrimalInfeasible())
		status = DSP_STAT_PRIM_INFEASIBLE;
	else if (si->isProvenDualInfeasible())
		status = DSP_STAT_DUAL_INFEASIBLE;
	else if (si->isPrimalObjectiveLimitReached())
		status = DSP_STAT_LIM_PRIM_OBJ;
	else if (si->isDualObjectiveLimitReached())
		status = DSP_STAT_LIM_DUAL_OBJ;
	else if (si->isIterationLimitReached())
		status = DSP_STAT_LIM_ITERorTIME;
	else if (si->isAbandoned())
		status = DSP_STAT_ABORT;
	else
		status = DSP_STAT_UNKNOWN;
}


#endif /* SRC_SOLVERINTERFACE_DSPOSI_H_ */