/*
 * DecBlkModel.h
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

#ifndef SRC_MODEL_DECBLKMODEL_H_
#define SRC_MODEL_DECBLKMODEL_H_

#include "Model/BlkModel.h"
#include "Model/DecModel.h"

class DecBlkModel: public DecModel {
public:

	/** default constructor */
	DecBlkModel();

	/** copy constructor */
	DecBlkModel(const DecBlkModel& rhs);

	/** default destructor */
	virtual ~DecBlkModel();

	/** get block model */
	BlkModel* blkPtr() {return blk_;}

public:

	std::vector<int> getCoupledSubproblemIndices(int j) {return blk_->getCoupledSubproblemIndices(j);}

	int getNumSubproblems() {return blk_->getNumBlocks() - 1;}

	int getNumCouplingRows() {return blk_->block(0)->getNumRows();}

	int getNumCouplingCols() {return blk_->block(0)->getNumCols();}

	int getNumSubproblemCouplingRows(int s) {return blk_->block(s+1)->getNumCouplingRows();}

	int getNumSubproblemCouplingCols(int s) {return blk_->block(s+1)->getNumCouplingCols();}

	const int * getSubproblemCouplingColIndices(int s) {return blk_->block(s+1)->getCouplingCols();}

	int getFullModelNumRows() {return blk_->getNumFullRows();}

	int getFullModelNumCols() {return blk_->getNumFullCols();}

	int getNumIntegers() {return blk_->getNumIntegers();}

	int getNumCouplingIntegers() {return blk_->block(0)->getNumIntegers();}

	//@{
	/** This function is used only for DD. */
	double evalLhsCouplingRow(int row, double ** solutions);
	double evalLhsCouplingRowSubprob(int row, int subprob, double * subprobSolution);
	//@}

	/**
	 * Retruns the coupling column objective coefficients
	 */
	const double * getCouplingColsObjs() {return blk_->block(0)->getObj();}

	double getCouplingRowLower(int row) {return blk_->block(0)->getRowLower()[row];}
	double getCouplingRowUpper(int row) {return blk_->block(0)->getRowUpper()[row];}

	char getSenseCouplingRow(int row);

	double getRhsCouplingRow(int row);

	bool nonanticipativity() {return false;}

	bool isStochastic() {return false;}
	bool isQCQP() {return false;}
	virtual void setDro(bool yes) { ; }
	bool isDro() {return false;}
	int getNumReferences() {return 0;}
	double getWassersteinSize() {return 0.0;}
	double getWassersteinDist(int i, int j) {return 0.0;}
	double getReferenceProbability(int i) {return 0.0;}

	DSP_RTN_CODE decompose(
		int size,                /**< [in] size of subproblem subset */
		int * subprobs,          /**< [in] subset of subproblems */
		int naux,                /**< [in] number of auxiliary columns */
		double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	DSP_RTN_CODE decompose(
		int size,                /**< [in] size of subproblem subset */
		int * subprobs,          /**< [in] subset of subproblems */
		int naux,                /**< [in] number of auxiliary columns */
		double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		CoinPackedMatrix *& qobj,/**< [out] quadratic objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	DSP_RTN_CODE decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */);

	/**
	 * This creates the subproblems with coupling columns; that is,
	 *   lb^k <= A^k x^k + B^k y^k <= ub^k
	 */
	DSP_RTN_CODE copySubprob(
		int subprob,             /**< [in] subproblem index */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix [A^k B^k] */
		double *& clbd,          /**< [out] column lower bounds of y */
		double *& cubd,          /**< [out] column upper bounds of y */
		char   *& ctype,         /**< [out] column types of y */
		double *& obj,           /**< [out] objective coefficients for y */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	DSP_RTN_CODE copyRecoProb(
		int scen,                      /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech,  /**< [out] technology matrix (A matrix) */
		CoinPackedMatrix *& mat_reco,  /**< [out] recourse matrix (B matrix) */
		double *& clbd_reco,           /**< [out] column lower bounds of y */
		double *& cubd_reco,           /**< [out] column upper bounds of y */
		char   *& ctype_reco,          /**< [out] column types of y */
		double *& obj_reco,            /**< [out] objective coefficients for y */
		double *& rlbd_reco,           /**< [out] row lower bounds */
		double *& rubd_reco,           /**< [out] row upper bounds */
		bool adjust_probability = true /**< not used */);

	DSP_RTN_CODE copyRecoProb(
		int scen,                     /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech, /**< [out] technology matrix (A matrix) */
		CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix (B matrix) */
		double *& clbd_reco,          /**< [out] column lower bounds of y */
		double *& cubd_reco,          /**< [out] column upper bounds of y */
		char   *& ctype_reco,         /**< [out] column types of y */
		double *& obj_reco,           /**< [out] objective coefficients for y */
		CoinPackedMatrix *& qobj_reco_coupling,/**< [out] coupling quadratric coefficients (y^2}*/
		CoinPackedMatrix *& qobj_reco_ncoupling, /**< [out] non-coupling quadratic coefficients (xy) */
		double *& rlbd_reco,          /**< [out] row lower bounds */
		double *& rubd_reco           /**< [out] row upper bounds */);

	DSP_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);

	DSP_RTN_CODE getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		CoinPackedMatrix *& qobj,/**< [out] quadratic objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */);
	void __printData();

protected:

	BlkModel* blk_; /**< block model pointer */
};

#endif /* SRC_MODEL_DECBLKMODEL_H_ */
