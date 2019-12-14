/**
 * test_de.cpp
 * 
 * 12/05/2019
 * Kibaek Kim
 */

#include "DspApiEnv.h"
#include "Model/DecTssModel.h"
#include "Solver/DualDecomp/DdDriverSerial.h"
#include "test_utils.hpp"

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        printf("Not enough arguments.\n");
        return 1;
    }
    else
    {
        DspApiEnv * env = new DspApiEnv;
        env->model_ = new DecTssModel;
        env->par_->setIntParam("DD/ITER_LIM", 3);

        TssModel * tss = dynamic_cast<TssModel*>(env->model_);
        tss->readSmps(argv[1]);
        
        env->solver_ = new DdDriverSerial(env->model_, env->par_, env->message_);
        env->solver_->init();
        dynamic_cast<DdDriverSerial*>(env->solver_)->run();
        env->solver_->finalize();

        double bestprimobj = env->solver_->getBestPrimalObjective();
        double bestdualobj = env->solver_->getBestDualObjective();
        double relgap = env->solver_->getRelDualityGap();
        printf("best primal bound: %.10f\n", bestprimobj);
        printf("best dual bound  : %.10f\n", bestdualobj);
        printf("relative gap (%%) : %f\n", relgap*100);

        FREE_PTR(env);

        double known = atof(argv[2]);
        printf("Known objective value: %+f\n", known);
        if (!isapprox(bestdualobj, known, bestprimobj)) {
            return 1;
        }
    }

    return 0;
}
