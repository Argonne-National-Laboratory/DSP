/*
 * SCIPconshdlrIntBendersWorker.cpp
 */

// #define DSP_DEBUG
#include "Utility/DspMessage.h"
#include "Solver/Benders/SCIPconshdlrIntBendersWorker.h"
#include "Solver/Benders/BdMW.h"

void SCIPconshdlrIntBendersWorker::setDecModel(DecModel *model)
{
	model_ = model;
	initialize(model_->getNumSubproblems());
	create_distsepa_problem();
}

SCIP_RETCODE SCIPconshdlrIntBendersWorker::evaluateRecourse(
	SCIP *scip,	   /**< [in] scip pointer */
	SCIP_SOL *sol, /**< [in] solution to evaluate */
	double *values /**< [out] evaluated recourse values */)
{
	/**< current solution */
	SCIP_Real *vals = NULL;

	/** recourse values collected from ranks but not sorted */
	double *recourse = NULL;

	/** Tell workers to evaluate recourse */
	int message = BdMW::MASTER_EVALUATES_RECOURSE;
	MPI_Bcast(&message, 1, MPI_INT, 0, comm_);

	/** allocate memory */
	SCIP_CALL(SCIPallocMemoryArray(scip, &vals, nvars_));

	/** get current solution */
	SCIP_CALL(SCIPgetSolVals(scip, sol, nvars_, vars_, vals));

	/** Send solutions to the workers */
	MPI_Bcast(vals, nvars_, MPI_DOUBLE, 0, comm_);

	/** check return value */
	int ret = 0, allret = 0;
	MPI_Allreduce(&ret, &allret, 1, MPI_INT, MPI_SUM, comm_);
	if (allret > 0)
		return SCIP_ERROR;

	/** receive recourse evaluations */
	recourse = new double[model_->getNumSubproblems()];
	MPI_Gatherv(NULL, 0, MPI_DOUBLE, recourse, recvcounts_, displs_, MPI_DOUBLE, 0, comm_);

	/** sort by index */
	int j = 0;
	for (int i = 1; i < comm_size_; ++i)
	{
		for (int s = i - 1; s < model_->getNumSubproblems(); s += comm_size_ - 1)
		{
			DSPdebugMessage("recourse[%d] = %e\n", s, recourse[j]);
			values[s] = recourse[j++];
		}
	}

	FREE_ARRAY_PTR(recourse)

	/** free memory */
	SCIPfreeMemoryArray(scip, &vals);

	return SCIP_OKAY;
}
