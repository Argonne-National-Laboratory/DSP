/**
 * test_de.cpp
 * 
 * 12/05/2019
 * Kibaek Kim
 */

#include "DspApiEnv.h"
#include "Utility/DspMpi.h"
#include "Model/DecTssModel.h"
#include "Solver/Benders/BdDriverMpi.h"
#include "test_utils.hpp"

const MPI_Comm comm = MPI_COMM_WORLD;

void distribute_procs(int nsubprobs, std::vector<int>& proc_idx_set) {
    int comm_rank, comm_size;
	MPI_Comm_rank(comm, &comm_rank);
	MPI_Comm_size(comm, &comm_size);

	// empty block IDs
	proc_idx_set.clear();

	// DSP is further parallelized with comm_size > nsubprobs.
	int modrank = comm_rank % nsubprobs;

    // assign sub-blocks in round-robin fashion
    for (int s = modrank; s < nsubprobs; s += comm_size) {
        proc_idx_set.push_back(s);
    }
}

int main(int argc, char* argv[])
{
    int ret = 0;

    if (argc < 3)
    {
        printf("Not enough arguments.\n");
        ret = 1;
    }
    else
    {
        int comm_rank, comm_size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(comm, &comm_rank);
        MPI_Comm_size(comm, &comm_size);

        DspApiEnv * env = new DspApiEnv;
        env->model_ = new DecTssModel;

        TssModel * tss = dynamic_cast<TssModel*>(env->model_);
        tss->readSmps(argv[1]);
        
        /** assign processes */
        std::vector<int> proc_idx_set;
        distribute_procs(env->model_->getNumSubproblems(), proc_idx_set);
        env->par_->setIntPtrParamSize("ARR_PROC_IDX", (int) proc_idx_set.size());
        for (int i = 0; i < proc_idx_set.size(); ++i)
            env->par_->setIntPtrParam("ARR_PROC_IDX", i, proc_idx_set[i]);
        
        env->solver_ = new BdDriverMpi(env->model_, env->par_, env->message_, comm);
        env->par_->setIntParam("BD/DD/ITER_LIM", 1);
        env->solver_->setTimeLimit(10.0);
        env->solver_->init();
        dynamic_cast<BdDriverMpi*>(env->solver_)->run();
        env->solver_->finalize();

        if (comm_rank == 0) {
            double bestprimobj = env->solver_->getBestPrimalObjective();
            double bestdualobj = env->solver_->getBestDualObjective();
            double relgap = env->solver_->getRelDualityGap();
            printf("best primal bound: %e\n", bestprimobj);
            printf("best dual bound  : %e\n", bestdualobj);
            printf("relative gap (%%) : %f\n", relgap*100);

            double known = atof(argv[2]);
            printf("Known objective value: %+f\n", known);
            if (!isapprox(bestdualobj, known, bestprimobj)) {
                ret = 1;
            }
        }
        MPI_Bcast(&ret, 1, MPI_INT, 0, comm);

        FREE_PTR(env);

        MPI_Finalize();
    }

    return ret;
}
