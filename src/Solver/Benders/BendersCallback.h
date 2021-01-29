#ifndef BENDERSCALLBACK_H_
#define BENDERSCALLBACK_H_

#include "Model/DecModel.h"
#include "Solver/Benders/BdSub.h"
#include "OsiCuts.hpp"
#include "SolverInterface/DspOsi.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiCpx.h"

struct cbdata{
	/** empty for now */
};

class BendersCallback{
    public: 
		/** default constructor */
		BendersCallback(DspOsi *osi, const char *name, int sepapriority);

		/** default constructor */
		virtual ~BendersCallback();

		/** set model pointer */
		virtual void setDecModel(DecModel * model) {model_ = model;}

        virtual void setBdSub(BdSub * bdsub);

		int static BendersCut(void *cbdata, int cbwhere);

		virtual DSP_RTN_CODE setOriginalVariables(
			int nvars,        /**< number of original variables, including auxiliary variables */
			int         naux  /**< number of auxiliary variables */);

		virtual DSP_RTN_CODE addBenderscut(DspOsi *osi);
        
    protected:

		virtual DSP_RTN_CODE generate_Benders(DspOsi * osi, OsiCuts *cs);
        /** generate Benders cuts */
        virtual void generateCuts(
		    int size,  /**< [in] size of x */
		    double *x, /**< [in] master solution */
		    OsiCuts *cuts /**< [out] cuts generated */);
    
        /** generate Benders cuts */
	    virtual void aggregateCuts(
			double ** cutvec, /**< [in] cut vector */
			double *  cutrhs, /**< [in] cut right-hand side */
			OsiCuts * cuts    /**< [out] cuts generated */);
	
	protected:
		DspOsi * osi_;					/** pointer to solver */
		DecModel *  model_;            /**< DecModel object */
		BdSub *     bdsub_;            /**< pointer to cut generator */
		int         nvars_;            /**< number of original variables */
		int         naux_;             /**< number of auxiliary variables */
		double*     probability_;      /**< array of probability */
};

#endif /** BENDERSCALLBACK_H_ */