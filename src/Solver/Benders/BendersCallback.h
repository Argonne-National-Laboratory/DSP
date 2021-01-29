#include "Model/DecModel.h"
#include "Solver/Benders/BdSub.h"
#include "OsiCuts.hpp"
#include "SolverInterface/DspOsi.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiCpx.h"

class BendersCallback:{
    public: 
		/** default constructor */
		BendersCallback(DspOsi *osi, const char *name, int sepapriority);

		/** default constructor */
		virtual ~BendersCallback();

        virtual void setBdSub(BdSub * bdsub);

		virtual int static Benderscut();

		virtual DSP_RTN_CODE setOriginalVariables(
			int nvars,        /**< number of original variables, including auxiliary variables */
			double   ** vars, /**< original variables, including auxiliary variables */
			int         naux  /**< number of auxiliary variables */);
        
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

		DecModel *  model_;            /**< DecModel object */
		BdSub *     bdsub_;            /**< pointer to cut generator */
		int         nvars_;            /**< number of original variables */
		double *	vars_;             /**< pointer array to original variables */
		int         naux_;             /**< number of auxiliary variables */
		double*     probability_;      /**< array of probability */
}