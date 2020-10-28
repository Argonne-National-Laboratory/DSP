/*
 * StoModel.h
 *
 *  Created on: Sep 22, 2014
 *      Author: kibaekkim
 */

#ifndef STOMODEL_H_
#define STOMODEL_H_

#include <map>
/** Coin */
#include "CoinTime.hpp"
#include "SmiScnModel.hpp"
/** Dsp */
#include "Utility/DspTypes.h"
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"


/*
 * This class is a wrapper for SMI, which reads and writes SMPS files.
 * A key justification of this class is to support various decomposition
 * methods by providing proper model structure.
 *
 * TODO: Better to derive SmiScnModel class?
 * TODO: Multi-stage stochastic optimization
 */
class StoModel {

	typedef std::map<int,int> StoScenMap;

public:

	/** default constructor */
	StoModel();

	/** copy constructor */
	StoModel(const StoModel & rhs);

	/** default destructor */
	virtual ~StoModel();

	/** read SMPS files */
	DSP_RTN_CODE readSmps(const char * filename);

	/** read DRO file */
	DSP_RTN_CODE readDro(const char * filename);

	void __printData();

public:

	/** get number of stages */
	int getNumStages() {return nstgs_;}

	/** get number of scenarios */
	int getNumScenarios() const {return nscen_;}

	/** get probability */
	const double * getProbability() const {return prob_;}

	/** get number of rows for a given stage */
	int getNumRows(int stage) const {return nrows_[stage];}

	/** get number of columns for a given stage */
	int getNumCols(int stage) const {return ncols_[stage];}

	/** get number of integer variables for a given stage */
	int getNumIntegers(int stage) const {return nints_[stage];}

	/** get number of integer variables in core */
	int getNumCoreIntegers() const {return nints_core_;}

	/** get objective function coefficients for a given stage */
	const double * getObjCore(int stage) {return obj_core_[stage];}

	const CoinPackedVector * getObjScenario(int scenario) {return obj_scen_[scenario];}

	/** get quadratic objective function coefficients for a given stage */
	const CoinPackedMatrix * getQuadraticObjCore(int stage) {return qobj_core_[stage];}

	const CoinPackedMatrix * getQuadraticObjScenario(int scenario) {return qobj_scen_[scenario];}

	/** get column type for a given stage */
	const char * getCtypeCore(int stage) {return ctype_core_[stage];}

	/** get initial solutions */
	const Solutions getInitialSolutions() {return init_solutions_;}

	/** set initial solutions */
	void setSolution(
			int      size,    /**< size of array */
			double * solution /**< solution */);

#if 0
	/** add branching object */
	void addBranchingHyperplane(int nzcnt, int * indices, double * values, int priority);
#endif

public:

	/** split core matrix row for a given stage */
	CoinPackedVector * splitCoreRowVec(
			int i,  /**< row index */
			int stg /**< stage */);

	/** copy core column lower bounds */
	void copyCoreColLower(double * clbd, int stg);

	/** copy core column upper bounds */
	void copyCoreColUpper(double * cubd, int stg);

	/** copy core objective coefficients */
	void copyCoreObjective(double * obj, int stg);

	/** copy core quadratic objective coefficients 
	 *  for coupling terms, do not shift indices, start from x,
	 *  for non-coupling terms, shift indices, start from y,
	 *  only for two-stage problem
	*/
	void copyCoreQuadrativeObjective(
		CoinPackedMatrix *&qobj_coupling, 
		CoinPackedMatrix *&qobj_ncoupling, 
		int stg);

	/** copy core column types */
	void copyCoreColType(char * ctype, int stg);

	/** copy core row lower bounds */
	void copyCoreRowLower(double * rlbd, int stg);

	/** copy core row upper bounds */
	void copyCoreRowUpper(double * rubd, int stg);

	/** combine random matrix row for a given scenario */
	void combineRandRowVec(
			CoinPackedVector * row, /**< core row vector */
			int i,                  /**< row index */
			int scen                /**< scenario index */);

	/** combine random matrix row for a given stage and scenario */
	void combineRandRowVec(
			CoinPackedVector * row, /**< core row vector */
			int i,                  /**< row index */
			int stg,                /**< stage index */
			int scen                /**< scenario index */);

	/** combine random column lower bounds */
	void combineRandColLower(double * clbd, int stg, int scen);

	/** combine random column upper bounds */
	void combineRandColUpper(double * cubd, int stg, int scen);

	/** combine random objective coefficients */
	void combineRandObjective(double * obj, int stg, int scen, bool adjustProbability = true);

	/** combine random quadratic objective coefficients */
	void combineRandQuadraticObjective(CoinPackedMatrix * &qobj_coupling, CoinPackedMatrix * &qobj_ncoupling, int stg, int scen, bool adjustProbability = true);

	/** combine random row lower bounds */
	void combineRandRowLower(double * rlbd, int stg, int scen);

	/** combine random row upper bounds */
	void combineRandRowUpper(double * rubd, int stg, int scen);

	/** shift vector indices by offset */
	void shiftVecIndices(
			CoinPackedVector * vec, /**< vector to shift indices */
			int offset,             /**< offset by which indices are shifted */
			int start = 0           /**< index only after which indices are shifted */);

	/** shift vector indices by offset */
	void shiftVecIndices(
			int size,     /**< size of vecind */
			int * vecind, /**< vector indices to shift */
			int offset,   /**< offset by which indices are shifted */
			int start = 0 /**< index only after which indices are shifted */);

	virtual bool isQCQP() {return isQCQP_;}

	// The following functions are for distributionally robust variant.
	// TODO: Better to create a new inhereted class?
	virtual void setDro(bool yes) { isdro_ = yes; }
	virtual bool isDro() {return isdro_;}
	virtual int getNumReferences() {return nrefs_;}
	virtual double getWassersteinSize() {return wass_eps_;}
	virtual double getWassersteinDist(int i, int j);
	virtual double getReferenceProbability(int i);

protected:

	/*
	 * Stage level data
	 */
	int nscen_;     /**< number of scenarios */
	int nstgs_;     /**< number of stages */
	int * nrows_;   /**< array of the number of rows for each stage */
	int * ncols_;   /**< array of the number of columns for each stage */
	int * nints_;   /**< number of integer variables for each stage */
	int * rstart_;  /**< array of row start indices with respect to core model */
	int * cstart_;  /**< array of column start indices with respect to core model */

	/*
	 * Core data.
	 */
	int nrows_core_;                /**< number of rows in core */
	int ncols_core_;                /**< number of columns in core */
	int nints_core_;                /**< number of integer variables in core */
	CoinPackedVector ** rows_core_; /**< rows in core matrix */
	double ** clbd_core_;           /**< column lower bounds for each stage */
	double ** cubd_core_;           /**< column upper bounds for each stage */
	double ** obj_core_;            /**< objective coefficients for each stage */
	CoinPackedMatrix ** qobj_core_; /**< quadratic objecitve coefficients for each stage */
	double ** rlbd_core_;           /**< row lower bounds for each stage */
	double ** rubd_core_;           /**< row upper bounds for each stage */
	char **   ctype_core_;          /**< column types for each stage */

	/*
	 * Random data only (no core data)
	 * TODO: This must assume two-stage programs.
	 */
	double * prob_;                 /**< array of scenario probability */
	CoinPackedMatrix ** mat_scen_;  /**< scenario matrix */
	CoinPackedVector ** clbd_scen_; /**< column lower bounds for each scenario */
	CoinPackedVector ** cubd_scen_; /**< column upper bounds for each scenario */
	CoinPackedVector ** obj_scen_;  /**< objective coefficients for each scenario */
	CoinPackedMatrix ** qobj_scen_; /**< quadratic objective coefficients for each scenario */
	CoinPackedVector ** rlbd_scen_; /**< row lower bounds for each scenario */
	CoinPackedVector ** rubd_scen_; /**< row upper bounds for each scenario */

	StoScenMap scen2stg_; /** map from scenario to stage */

	Solutions init_solutions_; /**< initial solutions */

	bool fromSMPS_; /**< problem was read from SMPS files? */

	bool isdro_;                /**< is this distributionally robust? */
	bool isQCQP_ = 0;			/**< quadratic terms in the problem? */
	int nrefs_;                 /**< number of reference scenarios for DRO */
	double wass_eps_;           /**< size of the Wasserstein ball */
	double ** wass_dist_;       /**< Wasserstein distances between two realizations */
	double * refs_probability_; /** probability vector of references */

public:

#if 0
	vector<pair<CoinPackedVector*,int> > branchingHyperplanes_; /**< branching hyper-planes */
#endif

};

#endif /* STOMODEL_H_ */
