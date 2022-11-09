/*
 * StoModel.h
 *
 *  Created on: Sep 22, 2014
 *      Author: kibaekkim
 */

#ifndef STOMODEL_H_
#define STOMODEL_H_

/** Coin */
#include "SmiScnModel.hpp"
#include "SmiScenarioTree.hpp"
#include "CoinTime.hpp"
/** Dsp */
#include "Utility/DspTypes.h"
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"

#include <vector>

struct StageData {
	int nrows_;   /**< array of the number of rows for each stage */
	int ncols_;   /**< array of the number of columns for each stage */
	int nints_;   /**< number of integer variables for each stage */
	int rstart_;  /**< array of row start indices with respect to core model */
	int cstart_;  /**< array of column start indices with respect to core model */
	double * clbd_core_;           /**< column lower bounds for each stage */
	double * cubd_core_;           /**< column upper bounds for each stage */
	double * obj_core_;            /**< objective coefficients for each stage */
	// CoinPackedMatrix * qobj_core_; /**< quadratic objecitve coefficients for each stage */
	double * rlbd_core_;           /**< row lower bounds for each stage */
	double * rubd_core_;           /**< row upper bounds for each stage */
	char *   ctype_core_;          /**< column types for each stage */
	CoinPackedVector ** rows_core_; /**< rows in core matrix */
	// QuadRowData * qc_row_core_;		/**< parameters for quadratic rows in core: current version only accept noncoupling quadratic rows */
	
	/** default constructor */
	StageData();

	/** copy constructor */
	StageData(const StageData & rhs);

	/** default destructor */
	virtual ~StageData();
};

class DspScnNode 
{
public:
	/** default constructor */
    DspScnNode();

	/** copy constructor */
	DspScnNode(const DspScnNode & rhs);

	/** default destructor */
	virtual ~DspScnNode();


	void addParent(DspScnNode* parent){ parent_ = parent;}
	void addChild(DspScnNode* child){ children_.push_back(child);}
	 /*
	 * Random data only (no core data)
	 */
	
	int stg_;						/** stage number*/
	double prob_;                 	/**< scenario probability */
	DspScnNode * parent_;
	std::vector<DspScnNode* > children_;
	CoinPackedMatrix * mat_scen_;  	/**< scenario matrix */
	double * clbd_scen_; 	/**< column lower bounds for each scenario */
	double * cubd_scen_; 	/**< column upper bounds for each scenario */
	double * obj_scen_;  	/**< objective coefficients for each scenario */
	// CoinPackedMatrix * qobj_scen_; /**< quadratic objective coefficients for each scenario */
	double * rlbd_scen_; 	/**< row lower bounds for each scenario */
	double * rubd_scen_; 	/**< row upper bounds for each scenario */
	// QuadRowData * qc_row_scen_;		/**< parameters for quadratic rows in scenarios: current version only accept noncoupling quadratic rows */
};

class StoModel {

public:

    /** default constructor */
	StoModel();

	/** copy constructor */
	StoModel(const StoModel & rhs);

	void copyAllSubNodes(const DspScnNode & rhs, DspScnNode* dsproot);

	/** default destructor */
	virtual ~StoModel();

    /** read SMPS files */
	DSP_RTN_CODE readSmps(const char * filename);

	/** read DRO file */
    // deal with dro later
	//DSP_RTN_CODE readDro(const char * filename);

	/** construct a map that maps variable names to their indices */
	bool mapVarnameIndex(map<string, int> &map_varName_index, const char * corefilename);

	/** read quadratic data file, extending the smps file */
    // deal with quad later
	//DSP_RTN_CODE readQuad(const char * smps, const char * filename);

	void __printData();

	/* print quadratic rows of scenario s, if s == -1, print quadratic rows in core */
    // deal with quad later
	//DSP_RTN_CODE printQuadRows (const int s);
	//DSP_RTN_CODE printQuadRows (const QuadRowData *qc);
public:

	/** get number of stages */
    int getNumStages() {return nstgs_;}

	/** get number of scenarios */
	int getNumScenarios() const {return nscen_;}

	/** get probability */
    const double getProbability(DspScnNode* node) {return node->prob_;}

	/** get number of rows for a given stage */
	int getNumRows(int stage) {return stage_data_[stage]->nrows_;}

	/** get number of columns for a given stage */
	int getNumCols(int stage) {return stage_data_[stage]->ncols_;}

	/** get number of integer variables for a given stage */
	int getNumIntegers(int stage) {return stage_data_[stage]->nints_;}

	/** get number of integer variables in core */
    // not used?
	int getNumCoreIntegers() {return nints_core_;}

	/** get number of quadratic constraints in core */
    // deal with quad later
	//int getNumCoreQRows() {return qc_row_core_->nqrows;}
	
	/** get number of quadratic constraints of a scenario*/
    // deal with quad later
	//int getNumScenQRows(int scen) {return qc_row_scen_[scen]->nqrows;}
	
	/** get objective function coefficients for a given stage */
	const double * getObjCore(int stage) {return stage_data_[stage]->obj_core_;}

	const double * getObjNode(DspScnNode* node);

	/** get quadratic objective function coefficients for a given stage */
    // deal with quad later
	//const CoinPackedMatrix * getQuadraticObjCore(int stage) {return qobj_core_[stage];}

	//const CoinPackedMatrix * getQuadraticObjScenario(int scenario) {return qobj_scen_[scenario];}

	/** get column type for a given stage */
	const char * getCtypeCore(int stage) {return stage_data_[stage]->ctype_core_;}

	/** get initial solutions */
    // not used?
	const Solutions getInitialSolutions() {return init_solutions_;}

	/** get core coefficeints for a given stage */
	const CoinPackedVector * getRowCore(int stage, int i) {return stage_data_[stage]->rows_core_[i];}

	/** get parameters for quadratic constraints in core*/
    // deal with quad later
	//QuadRowData * getQuaraticsRowCore() const {return qc_row_core_;}

	/** get parameters for quadratic constraints in a scenario */
    // deal with quad later
	//QuadRowData * getQuaraticsRowScenario(int s) const {return qc_row_scen_[s];}

	//bool hasQuadraticRowCore() const {return qc_row_core_ != NULL ? true : false;};
	//bool hasQuadraticRowScenario() const {return qc_row_scen_ != NULL ? true : false;};
	//bool hasQuadraticRow() const {return (hasQuadraticRowCore() || hasQuadraticRowScenario());}

	/** set probability */
	void setProbability(DspScnNode* node, double probability);

	/** set initial solutions */
	void setSolution(
			int      size,    /**< size of array */
			double * solution /**< solution */);

	/** 
	 * Set the Wasserstein ambiguity set for distributionally robust optimization.
	 * This should be used for stochastic programming models, where the probabilities
	 * are used as the empirical references of the Wasserstein distance.
	 * 
	 * @param lp_norm use the distance norm of p
	 * @param eps maximum distance size
	 * @return DSP_RTN_OK if no error
	 */

    // deal with dro later
	//DSP_RTN_CODE setWassersteinAmbiguitySet(double lp_norm, double eps);

	/** 
	 * Nomalize probability vector
	 */
    // deal with dro later
	//void normalizeProbability();

#if 0
	/** add branching object */
	void addBranchingHyperplane(int nzcnt, int * indices, double * values, int priority);
#endif

private:
	int counter; /**< for verification of core matrix order*/

	/** load stage from core */
	void addStage(
		SmiCoreData * core,
		int s
	);

	DspScnNode* createNode(
		SmiCoreData * core,
		SmiScnNode * scnnode
	);

	void addAllSubNodes(
		SmiCoreData * core,
		SmiTreeNode<SmiScnNode* > * root,
		DspScnNode* dsproot
	);

	/** split given row vec*/
	CoinPackedVector * splitRowVec(
		SmiNodeData * node,
		SmiCoreData * core,
		int i /**< row index */
	);
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
	// void copyCoreQuadraticObjective(
	// 	CoinPackedMatrix *&qobj_coupling,
	// 	CoinPackedMatrix *&qobj_ncoupling,
	// 	int stg);

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
			DspScnNode* node        /**< node index */);

	/** combine random matrix row for a given stage and scenario */
	void combineRandRowVec(
			CoinPackedVector * row, /**< core row vector */
			int i,                  /**< row index */
			int stg,                /**< stage index */
			DspScnNode* node        /**< node index */);

	/** combine random column lower bounds */
	void combineRandColLower(double * clbd, DspScnNode* node);

	/** combine random column upper bounds */
	void combineRandColUpper(double * cubd, DspScnNode* node);

	/** combine random objective coefficients */
	void combineRandObjective(double * obj, DspScnNode* node, double adjustProbability = 1.0);

	/** combine random quadratic objective coefficients */
	// void combineRandQuadraticObjective(CoinPackedMatrix * &qobj_coupling, CoinPackedMatrix * &qobj_ncoupling, int stg, int scen, bool adjustProbability = true);

	/** combine random row lower bounds */
	void combineRandRowLower(double * rlbd, DspScnNode* node);

	/** combine random row upper bounds */
	void combineRandRowUpper(double * rubd, DspScnNode* node);

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

	// The following functions are for distributionally robust variant.
	// TODO: Better to create a new inhereted class?
	// virtual void setDro(bool yes) { isdro_ = yes; }
	// virtual bool isDro() {return isdro_;}
	// virtual int getNumReferences() {return nrefs_;}
	// virtual double getWassersteinSize() {return wass_eps_;}
	// virtual double getWassersteinDist(int i, int j);
	// virtual double getReferenceProbability(int i);

protected:

	// /*
	//  * Stage level data
	//  */
	int nstgs_;     /**< number of stages */
	int nnodes_;	/**< number of nodes */
	int nscen_;     /**< number of scenarios */

	// /*
	//  * Core data.
	//  */
	int nrows_core_;                /**< number of rows in core */
	int ncols_core_;                /**< number of columns in core */
	int nints_core_;                /**< number of integer variables in core */

	StageData** stage_data_;
    std::vector<DspScnNode* > node_data_;

	Solutions init_solutions_; /**< initial solutions */

	bool fromSMPS_; /**< problem was read from SMPS files? */

	// bool isdro_;				/**< is this distributionally robust? */
	// int nrefs_;                 /**< number of reference scenarios for DRO */
	// double wass_eps_;           /**< size of the Wasserstein ball */
	// double ** wass_dist_;       /**< Wasserstein distances between two realizations */
	// double * refs_probability_; /** probability vector of references */

public:

#if 0
	vector<pair<CoinPackedVector*,int> > branchingHyperplanes_; /**< branching hyper-planes */
#endif

};



#endif /** STOMODEL_H_ */