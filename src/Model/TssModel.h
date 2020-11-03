/*
 * TssModel.h
 *
 *  Created on: Sep 23, 2014
 *      Author: kibaekkim
 */

#ifndef TSSMODEL_H_
#define TSSMODEL_H_

#define TSSMODEL_DEBUG

#include "StoModel.h"
#include "DetModel.h"

class TssModel: public StoModel {
public:

	/** default constructor */
	TssModel();

	/** copy constructor */
	TssModel(const TssModel & rhs);

	/** default destructor */
	virtual ~TssModel();

public:

	/**
	 * This creates recourse objective function for a given scenario index.
	 */
	DSP_RTN_CODE copyRecoObj(
			int scen,                     /**< [in] scenario index */
			double *& obj_reco,           /**< [out] objective coefficients */
			bool adjustProbability = true);

	DSP_RTN_CODE copyRecoObj(
			int scen,                     /**< [in] scenario index */
			double *& obj_reco,           /**< [out] objective coefficients */
			CoinPackedMatrix *& qobj_reco_coupling,/**< [out] coupling quadratric coefficients (y^2}*/
			CoinPackedMatrix *& qobj_reco_ncoupling, /**< [out] non-coupling quadratic coefficients (xy) */
			bool adjustProbability = true);
	/** for C API functions */

public:

	/** set number of scenarios */
	DSP_RTN_CODE setNumberOfScenarios(int nscen);

	/** set dimensions */
	DSP_RTN_CODE setDimensions(
			const int ncols1, /**< number of first-stage columns */
			const int nrows1, /**< number of first-stage rows */
			const int ncols2, /**< number of second-stage columns */
			const int nrows2  /**< number of second-stage rows */);

	/** load first-stage problem */
	DSP_RTN_CODE loadFirstStage(
			const CoinBigIndex * start, /**< start index for each row */
			const int *          index, /**< column indices */
			const double *       value, /**< constraint elements */
			const double *       clbd,  /**< column lower bounds */
			const double *       cubd,  /**< column upper bounds */
			const char *         ctype, /**< column types */
			const double *       obj,   /**< objective coefficients */
			const double *       rlbd,  /**< row lower bounds */
			const double *       rubd   /**< row upper bounds */);

	/** load first-stage problem with quadratic objectives */
	DSP_RTN_CODE loadFirstStage(
			const CoinBigIndex * start, /**< start index for each row */
			const int *          index, /**< column indices */
			const double *       value, /**< constraint elements */
			const double *       clbd,  /**< column lower bounds */
			const double *       cubd,  /**< column upper bounds */
			const char *         ctype, /**< column types */
			const double *       obj,   /**< objective coefficients */
			const int * 		 qobjrowindex, /**< quadratic objective row indices */
			const int *			 qobjcolindex, /**< quadratic objective column indices */
			const double *		 qobjvalue, /**< quadratic objective constraint elements value */
			const int 			 numq,  /**< number of quadratic terms */
			const double *       rlbd,  /**< row lower bounds */
			const double *       rubd   /**< row upper bounds */);

	/** load second-stage problem */
	DSP_RTN_CODE loadSecondStage(
			const int            s,     /**< scenario index */
			const double         prob,  /**< probability */
			const CoinBigIndex * start, /**< start index for each row */
			const int *          index, /**< column indices */
			const double *       value, /**< constraint elements */
			const double *       clbd,  /**< column lower bounds */
			const double *       cubd,  /**< column upper bounds */
			const char *         ctype, /**< column types */
			const double *       obj,   /**< objective coefficients */
			const double *       rlbd,  /**< row lower bounds */
			const double *       rubd   /**< row upper bounds */);

	DSP_RTN_CODE loadSecondStage(
			const int            s,     /**< scenario index */
			const double         prob,  /**< probability */
			const CoinBigIndex * start, /**< start index for each row */
			const int *          index, /**< column indices */
			const double *       value, /**< constraint elements */
			const double *       clbd,  /**< column lower bounds */
			const double *       cubd,  /**< column upper bounds */
			const char *         ctype, /**< column types */
			const double *       obj,   /**< objective coefficients */
			const int * 		 qobjrowindex, /**< quadratic objective row indices */
			const int *			 qobjcolindex, /**< quadratic objective column indices */
			const double *		 qobjvalue, /**< quadratic objective constraint elements value */
			const int 			 numq,  /**< number of quadratic terms */
			const double *       rlbd,  /**< row lower bounds */
			const double *       rubd   /**< row upper bounds */);
};

#endif /* TSSMODEL_H_ */
