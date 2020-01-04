/*
 * DwWorkerMpi.h
 *
 *  Created on: Nov 18, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_SOLVER_DANTZIGWOLFE_DWWORKERMPI_H_
#define SRC_SOLVER_DANTZIGWOLFE_DWWORKERMPI_H_

#include "Utility/DspMpi.h"
#include "Solver/DantzigWolfe/DwWorker.h"

class DwWorkerMpi: public DwWorker {
public:

	/** default constructor */
	DwWorkerMpi(
			DecModel * model,
			DspParams * par,
			DspMessage * message,
			MPI_Comm comm);

	/** default destructor */
	virtual ~DwWorkerMpi();

	/** generate variables */
	virtual DSP_RTN_CODE generateCols(
			int phase,                           /**< [in] phase of the master */
			const double* piA,                   /**< [in] piA */
			std::vector<int>& indices,           /**< [out] subproblem indices */
			std::vector<int>& statuses,          /**< [out] solution status */
			std::vector<double>& cxs,            /**< [out] solution times original objective coefficients */
			std::vector<double>& objs,           /**< [out] subproblem objective values */
			std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */);

	/** generate variables by fixing some of the varialbes (e.g., upper bounding for SMIP) */
	virtual DSP_RTN_CODE generateColsByFix(
			const double* x,                     /**< [in] solution to fix */
			std::vector<int>& indices,           /**< [out] subproblem indices */
			std::vector<int>& statuses,          /**< [out] solution status */
			std::vector<double>& objs,           /**< [out] subproblem objective values */
			std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */);

	/** set column bounds */
	virtual void setColBounds(int size, const int* indices, const double* lbs, const double* ubs);

	/** add row (for branching disjunction) */
	virtual void addRow(const CoinPackedVector* vec, double lb, double ub);

	/** reset time increment */
	virtual void resetTimeIncrement() {
		resetTimeIncrement_ = 1;
	}

	/** In this function, non-root processes receive signals
	 * from the root process and do proper processing. */
	virtual DSP_RTN_CODE receiver();

	enum {
		sig_generateCols = 0,
		sig_generateColsByFix,
		sig_setColBounds,
		sig_addRow,
		sig_terminate
	};

protected:

	MPI_Comm comm_;
	int comm_rank_;
	int comm_size_;

	int npiA_; /**< size of piA vector */
	double* piA_; /**< local piA vector */
	int resetTimeIncrement_;
};

#endif /* SRC_SOLVER_DANTZIGWOLFE_DWWORKERMPI_H_ */
