/*
 * DetModel.h
 *
 *  Created on: July 1, 2015
 *      Author: ctjandra
 */

#ifndef DETMODEL_H_
#define DETMODEL_H_

#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"

/**
 * Deterministic model. To be used with general decomposition models (DecModel) that are not based on stochastic
 * models.
 */
class DetModel {

public:

	/** default constructor */
	DetModel();

	DetModel(
		const CoinPackedMatrix * mat,   /**< constraint matrix */
		const double * clbd,            /**< column lower bounds */
		const double * cubd,            /**< column upper bounds */
		const char   * ctype,           /**< column types */
		const double * obj,             /**< objective coefficients */
		const double * rlbd,            /**< row lower bounds */
		const double * rubd             /**< row upper bounds */);

	DetModel(
		const CoinPackedMatrix * mat,   /**< constraint matrix */
		const double * clbd,            /**< column lower bounds */
		const double * cubd,            /**< column upper bounds */
		const char   * ctype,           /**< column types */
		const double * obj,             /**< objective coefficients */
		const CoinPackedMatrix * qobj,  /**< quadratic objective coefficients */
		const double * rlbd,            /**< row lower bounds */
		const double * rubd             /**< row upper bounds */);

	DetModel(
		const CoinBigIndex * start,     /**< start index for each row */
		const int    * index,           /**< column indices */
		const double * value,           /**< constraint elements */
		const int      numels,          /**< number of elements in index and value */
		const int      ncols,           /**< number of columns */
		const int      nrows,           /**< number of rows */
		const double * clbd,            /**< column lower bounds */
		const double * cubd,            /**< column upper bounds */
		const char   * ctype,           /**< column types */
		const double * obj,             /**< objective coefficients */
		const double * rlbd,            /**< row lower bounds */
		const double * rubd             /**< row upper bounds */);
	
	DetModel(
		const CoinBigIndex * start,     /**< start index for each row */
		const int    * index,           /**< column indices */
		const double * value,           /**< constraint elements */
		const int      numels,          /**< number of elements in index and value */
		const int      ncols,           /**< number of columns */
		const int      nrows,           /**< number of rows */
		const double * clbd,            /**< column lower bounds */
		const double * cubd,            /**< column upper bounds */
		const char   * ctype,           /**< column types */
		const double * obj,             /**< objective coefficients */
		const CoinBigIndex * qobjstart, /**< quadratic objective start index for each row */
		const int	 * qobjindex,		/**< quadratic objective index */
		const double * qobjvalue,		/**< quadratic objective coefficient*/
		const int	   qobjnumels,		/**< number of elements in qobj index and value */
		const double * rlbd,            /**< row lower bounds */
		const double * rubd             /**< row upper bounds */);

	/** copy constructor */
	DetModel(const DetModel & rhs);

	/** default destructor */
	virtual ~DetModel();

	/** create model */
	void createModel(
			const CoinPackedMatrix * mat, /**< constraint matrix */
			const double * clbd,          /**< column lower bounds */
			const double * cubd,          /**< column upper bounds */
			const char   * ctype,         /**< column types */
			const double * obj,           /**< objective coefficients */
			const double * rlbd,          /**< row lower bounds */
			const double * rubd           /**< row upper bounds */);

	/** create model */
	void createModel(
			const CoinPackedMatrix * mat,   /**< constraint matrix */
			const double * clbd,            /**< column lower bounds */
			const double * cubd,            /**< column upper bounds */
			const char   * ctype,           /**< column types */
			const double * obj,             /**< objective coefficients */
			const CoinPackedMatrix * qobj,  /**< quadratic objective coefficients */
			const double * rlbd,            /**< row lower bounds */
			const double * rubd             /**< row upper bounds */);

	/** create model */
	void createModel(
			const CoinBigIndex * start, /**< start index for each row */
			const int    * index,       /**< column indices */
			const double * value,       /**< constraint elements */
			const int      numels,      /**< number of elements in index and value */
			const int      ncols,       /**< number of columns */
			const int      nrows,       /**< number of rows */
			const double * clbd,        /**< column lower bounds */
			const double * cubd,        /**< column upper bounds */
			const char   * ctype,       /**< column types */
			const double * obj,         /**< objective coefficients */
			const double * rlbd,        /**< row lower bounds */
			const double * rubd         /**< row upper bounds */);
	
	/** create model */
	void createModel(
			const CoinBigIndex * start, /**< start index for each row */
			const int    * index,       /**< column indices */
			const double * value,       /**< constraint elements */
			const int      numels,      /**< number of elements in index and value */
			const int      ncols,       /**< number of columns */
			const int      nrows,       /**< number of rows */
			const double * clbd,        /**< column lower bounds */
			const double * cubd,        /**< column upper bounds */
			const char   * ctype,       /**< column types */
			const double * obj,         /**< objective coefficients */
			const CoinBigIndex * qobjstart, /**< quadratic objective start index for each row */
			const int	 * qobjindex,	/**< quadratic objective index */
			const double * qobjvalue,	/**< quadratic objective coefficient*/
			const int	   qobjnumels,  /**< number of elements in qobj index and value */
			const double * rlbd,        /**< row lower bounds */
			const double * rubd         /**< row upper bounds */);
	void __printData();

public:

	const CoinPackedMatrix * getConstraintMatrix() {return mat_;}

	/** get number of rows */
	int getNumRows() const {return mat_->getNumRows();}

	/** get number of columns */
	int getNumCols() const {return mat_->getNumCols();}

	/** get number of integer variables */
	int getNumIntegers() const {return nints_;}

	/** get column lower bounds */
	const double * getColLower() {return clbd_;}

	/** get column upper bounds */
	const double * getColUpper() {return cubd_;}

	/** get column types */
	const char * getCtype() {return ctype_;}

	/** get objective function coefficients */
	const double * getObj() {return obj_;}

	/** get full quadratic objective matrix */
	const CoinPackedMatrix * getQuadraticObjectiveMatrix() {return qobj_;}

	/** get row lower bounds */
	const double * getRowLower() {return rlbd_;}

	/** get row upper bounds */
	const double * getRowUpper() {return rubd_;}

protected:

	/** Model data */
	CoinPackedMatrix * mat_;   /**< constraint matrix */
	double * clbd_;            /**< column lower bounds */
	double * cubd_;            /**< column upper bounds */
	char   * ctype_;           /**< column types */
	double * obj_;             /**< objective coefficients */
	CoinPackedMatrix * qobj_;  /**< quadratic objective coefficients */
	double * rlbd_;            /**< row lower bounds */
	double * rubd_;            /**< row upper bounds */
	int nints_;                /**< number of integer variables */

};

#endif /* DETMODEL_H_ */
