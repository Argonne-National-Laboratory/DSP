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

/** quadratic row infomration for a scenario */
struct QuadRowData {
	int				nqrows; 		/** number of quadratic rows for a scenario  */
    int *         	linnzcnt;  		/** number of nonzero coefficients in the linear part of each constraint  */
    int *        	quadnzcnt;  	/** number of nonzero coefficients in the quadratic part of each constraint  */
	double *		rhs; 			/** constraint rhs of each constraint */
	int *			sense; 			/** constraint sense of each constraint */
	int **        	linind; 		/** indices for the linear part */
	double **      	linval; 		/** nonzero coefficient of the linear part */
	int **      	quadrow;  		/** indices for the quadratic part */
	int **      	quadcol;  		/** indices for the quadratic part */
	double **     	quadval; 		/** nonzero coefficient of the quadratic part */ 

	QuadRowData(){
		nqrows = 0;
		linnzcnt = NULL;
		quadnzcnt = NULL;
		rhs = NULL;
		sense = NULL;
		linind = NULL;
		linval = NULL;
		quadrow = NULL;
		quadcol = NULL;
		quadval = NULL;
	}
	QuadRowData(int nqrows_, int rstart, vector<char> sense_, vector<double> rhs_, vector<vector<int>> linind_, vector<vector<double>> linval_, vector<vector<int>> quadrow_, vector<vector<int>> quadcol_, vector<vector<double>> quadval_)
	{
		nqrows = nqrows_;
		if (nqrows > 0) 
		{
			linnzcnt = new int [nqrows];
			quadnzcnt = new int [nqrows];
			rhs = new double [nqrows];
			sense = new int [nqrows];
			linind = new int * [nqrows];
			linval = new double * [nqrows];
			quadrow = new int * [nqrows];
			quadcol = new int * [nqrows];
			quadval = new double * [nqrows];

			for (int i = 0; i < nqrows; i++) 
			{
				sense[i] = sense_[i+rstart];
				rhs[i] = rhs_[i+rstart];
				
				linnzcnt[i] = linval_[i+rstart].size();
				linind[i] = new int [linnzcnt[i]];
				linval[i] = new double [linnzcnt[i]];

				quadnzcnt[i] = quadrow_[i+rstart].size();
				quadrow[i] = new int [quadnzcnt[i]];
				quadcol[i] = new int [quadnzcnt[i]];
				quadval[i] = new double [quadnzcnt[i]];

				for (int j = 0; j < linnzcnt[i]; j++) 
				{
					linind[i][j] = linind_[i+rstart][j];
					linval[i][j] = linval_[i+rstart][j];
				}

				for (int j = 0; j < quadnzcnt[i]; j++) 
				{
					quadrow[i][j] = quadrow_[i+rstart][j];
					quadcol[i][j] = quadcol_[i+rstart][j];
					quadval[i][j] = quadval_[i+rstart][j];
				}
			}
		}
	}
	~QuadRowData(){
		if (nqrows > 0) {
			FREE_ARRAY_PTR(linnzcnt);
			FREE_ARRAY_PTR(quadnzcnt);
			FREE_ARRAY_PTR(rhs);
			FREE_ARRAY_PTR(sense);

			FREE_2D_ARRAY_PTR(nqrows, linind);
			FREE_2D_ARRAY_PTR(nqrows, linval);
			FREE_2D_ARRAY_PTR(nqrows, quadrow);
			FREE_2D_ARRAY_PTR(nqrows, quadcol);
			FREE_2D_ARRAY_PTR(nqrows, quadval);
		}
	}
};

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

	/** construct a map that maps variable names to their indices */
	bool mapVarnameIndex(map<string, int> &map_varName_index, const char * corefilename);

	/** read quadratic data file, extending the smps file */
	DSP_RTN_CODE readQuad(const char * smps, const char * filename, bool chg_to_socp = true);

	void __printData();

	/* print quadratic rows of scenario s, if s == -1, print quadratic rows in core */
	DSP_RTN_CODE printQuadRows (const int s);
	DSP_RTN_CODE printQuadRows (const QuadRowData *qc);

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

	/** get number of quadratic constraints in core */
	int getNumCoreQRows() {return qc_row_core_->nqrows;}
	
	/** get number of quadratic constraints of a scenario*/
	int getNumScenQRows(int scen) {return qc_row_scen_[scen]->nqrows;}
	
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

	/** get core coefficeints for a given stage */
	const CoinPackedVector * getRowCore(int i) {return rows_core_[i];}

	/** get parameters for quadratic constraints in core*/
	QuadRowData * getQuaraticsRowCore() const {return qc_row_core_;}

	/** get parameters for quadratic constraints in a scenario */
	QuadRowData * getQuaraticsRowScenario(int s) const {return qc_row_scen_[s];}

	DSP_RTN_CODE chgToSocp(vector<int> &qc_rstart);
	DSP_RTN_CODE getL(double * &Q, int quadnzcnt, int *quadcol, int *quadrow, double *quadval, vector<int> &indices, int &n);

	bool hasQuadraticRowCore() const {return qc_row_core_ != NULL ? true : false;};
	bool hasQuadraticRowScenario() const {return qc_row_scen_ != NULL ? true : false;};
	bool hasQuadraticRow() const {return (hasQuadraticRowCore() || hasQuadraticRowScenario());}

	/** set probability */
	void setProbability(double *probability);

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
	DSP_RTN_CODE setWassersteinAmbiguitySet(double lp_norm, double eps);

	/** 
	 * Nomalize probability vector
	 */
	void normalizeProbability();

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
	void copyCoreQuadraticObjective(
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

	// The following functions are for distributionally robust variant.
	// TODO: Better to create a new inhereted class?
	virtual void setDro(bool yes) { isdro_ = yes; }
	virtual bool isDro() {return isdro_;}
	virtual int getNumReferences() {return nrefs_;}
	virtual double getWassersteinSize() {return wass_eps_;}
	virtual double getWassersteinDist(int i, int j);
	virtual double getReferenceProbability(int i);

	virtual bool isQcp() {return qc_row_scen_ == NULL ? false : true;}

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
	QuadRowData * qc_row_core_;		/**< parameters for quadratic rows in core: current version only accept noncoupling quadratic rows */

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
	QuadRowData ** qc_row_scen_;		/**< parameters for quadratic rows in scenarios: current version only accept noncoupling quadratic rows */

	StoScenMap scen2stg_; /** map from scenario to stage */

	Solutions init_solutions_; /**< initial solutions */

	bool fromSMPS_; /**< problem was read from SMPS files? */

	bool isdro_;				/**< is this distributionally robust? */
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
