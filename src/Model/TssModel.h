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
	 * This routine decomposes the model based on inputs. If size = 0, then
	 * this results in a standard Benders decomposition structure. If size = 1, then
	 * this results in a dual decomposition structure.
	 */
	DSP_RTN_CODE decompose(
			int size,                /**< [in] size of scenario subset */
			int * scen,              /**< [in] subset of scenarios */
			int naux,                /**< [in] number of auxiliary columns */
			double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
			double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
			double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
			//CoinPackedMatrix * qobj_aux, /**< [in] quadratic objective coefficients for auxiliary columns */
			CoinPackedMatrix *& mat, /**< [out] constraint matrix */
			double *& clbd,          /**< [out] column lower bounds */
			double *& cubd,          /**< [out] column upper bounds */
			char   *& ctype,         /**< [out] column types */
			double *& obj,           /**< [out] objective coefficients */
			//CoinPackedMatrix * qobj, /**< [out] quadratic objective coefficients */
			double *& rlbd,          /**< [out] row lower bounds */
			double *& rubd           /**< [out] row upper bounds */);

	DSP_RTN_CODE decompose(
			int size,                /**< [in] size of scenario subset */
			int * scen,              /**< [in] subset of scenarios */
			int naux,                /**< [in] number of auxiliary columns */
			double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
			double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
			double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
			//CoinPackedMatrix * qobj_aux, /**< [in] quadratic objective coefficients for auxiliary columns */
			CoinPackedMatrix *& mat, /**< [out] constraint matrix */
			double *& clbd,          /**< [out] column lower bounds */
			double *& cubd,          /**< [out] column upper bounds */
			char   *& ctype,         /**< [out] column types */
			double *& obj,           /**< [out] objective coefficients */
			CoinPackedMatrix * &qobj,/**< [out] quadratic objective coefficients */
			double *& rlbd,          /**< [out] row lower bounds */
			double *& rubd           /**< [out] row upper bounds */);

	/**
	 * This constructs a deterministic equivalent form.
	 */
	DSP_RTN_CODE copyDeterministicEquivalent(
			CoinPackedMatrix *& mat, /**< [out] constraint matrix */
			double *& clbd,          /**< [out] column lower bounds */
			double *& cubd,          /**< [out] column upper bounds */
			char   *& ctype,         /**< [out] column types */
			double *& obj,           /**< [out] objective coefficients */
			CoinPackedMatrix *& qobj,/**< [out] quadratic objective coefficients */
			double *& rlbd,          /**< [out] row lower bounds */
			double *& rubd           /**< [out] row upper bounds */);
			
	/**
	 * Create a DetModel representing the deterministic equivalent of this model.
	 */
	DSP_RTN_CODE copyDeterministicEquivalent(
			DetModel *& det /**< [out] deterministic equivalent model */);

	/**
	 * This creates recourse problem structure for a given scenario index.
	 */
	DSP_RTN_CODE copyRecoProb(
			int scen,                     /**< [in] scenario index */
			CoinPackedMatrix *& mat_tech, /**< [out] technology matrix */
			CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix */
			double *& clbd_reco,          /**< [out] column lower bounds */
			double *& cubd_reco,          /**< [out] column upper bounds */
			char   *& ctype_reco,         /**< [out] column types */
			double *& obj_reco,           /**< [out] objective coefficients */
			//CoinPackedMatrix *& qobj_reco,/**< [out] quadratic objective coefficients */
			double *& rlbd_reco,          /**< [out] row lower bounds */
			double *& rubd_reco           /**< [out] row upper bounds */);

	DSP_RTN_CODE copyRecoProb(
			int scen,                     /**< [in] scenario index */
			CoinPackedMatrix *& mat_tech, /**< [out] technology matrix */
			CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix */
			double *& clbd_reco,          /**< [out] column lower bounds */
			double *& cubd_reco,          /**< [out] column upper bounds */
			char   *& ctype_reco,         /**< [out] column types */
			double *& obj_reco,           /**< [out] objective coefficients */
			CoinPackedMatrix *& qobj_reco,/**< [out] quadratic objective coefficients */
			double *& rlbd_reco,          /**< [out] row lower bounds */
			double *& rubd_reco           /**< [out] row upper bounds */);

	/**
	 * This creates recourse objective function for a given scenario index.
	 */
	DSP_RTN_CODE copyRecoObj(
			int scen,                     /**< [in] scenario index */
			double *& obj_reco,           /**< [out] objective coefficients */
			bool adjustProbability = true);

	DSP_RTN_CODE copyRecoQuadraticObj(
			int scen,                     /**< [in] scenario index */
			CoinPackedMatrix *& qobj_reco,/**< [out] quadratic objective coefficients */
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

	/** load first-stage problem */
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
