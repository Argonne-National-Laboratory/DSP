/*
 * QCModel.h
 *
 *  Created on: Oct 22, 2020
 *      Author: geunyeongbyeon
 */

#ifndef QCMODEL_H_
#define QCMODEL_H_

//#define QCMODEL_DEBUG

#include <vector>
#include <algorithm>
#include <fstream>
// #include <map>
/** Coin */
// #include "CoinTime.hpp"
// #include "CoinPackedMatrix.hpp"
// #include "CoinHelperFunctions.hpp"
/** Dsp */
#include "Utility/DspTypes.h"
#include "Utility/DspMacros.h"
#include "Utility/DspRtnCodes.h"

/** quadratic row infomration for a scenario */
struct QcRowDataScen {
	int			nqrows_; 			/** number of quadratic rows for a scenario  */
    int *         linnzcnt_;  		/** number of nonzero coefficients in the linear part of each constraint  */
    int *        	quadnzcnt_;  	/** number of nonzero coefficients in the quadratic part of each constraint  */
	double *		rhs_; 			/** constraint rhs of each constraint */
	int *			sense_; 		/** constraint sense of each constraint */
	int **        linind_; 			/** indices for the linear part */
	double **      linval_; 		/** nonzero coefficient of the linear part */
	int **      	quadrow_;  		/** indices for the quadratic part */
	int **      	quadcol_;  		/** indices for the quadratic part */
	double **     quadval_; 		/** nonzero coefficient of the quadratic part */ 
};

class QcModel {
public:

	/** default constructor */
	QcModel();

	/** copy constructor */
	QcModel(const QcModel & rhs);

	/** default destructor */
	virtual ~QcModel();

public:

	/** construct a map that maps variable names to their indices */
	bool mapVarnameIndex(map<string, int> &map_varName_index, const char * corefilename);
	
	/** read quadratic data file */
	DSP_RTN_CODE readQuad(const char * smps, const char * filename, int ncols_2nd_);

	/** set dimensions for quadratic constraints */
	DSP_RTN_CODE setQuadDimensions(int nscen);
	DSP_RTN_CODE setQuadDimensions(int s, int nqrows);
	DSP_RTN_CODE setQuadDimensions(const int nscen, int * nqrows);

	/** set file name */
	void setFileName(const char * smps) {filename_ = smps;} 

	/** add quadratic constraints to the second stage problem */
	DSP_RTN_CODE loadQuadraticConstrs (
			const int           s,     		/**< scenario index */
		const int 			nqrows,
        const int *         linnzcnt,  	/**< number of nonzero coefficients in the linear part of each constraint  */
        const int *        	quadnzcnt,  /**< number of nonzero coefficients in the quadratic part of each constraint  */
		const double *		rhs, 		/**< constraint rhs of each constraint */
		const int *			sense, 		/**< constraint sense of each constraint */
		const int *         linstart,  	/**< number of nonzero coefficients in the linear part of each constraint  */
		const int *         linind, 	/**< indices for the linear part */
		const double *      linval, 	/**< nonzero coefficient of the linear part */
		const int *        	quadstart,  /**< number of nonzero coefficients in the quadratic part of each constraint  */
		const int *       	quadrow,  	/**< indices for the quadratic part */
		const int *       	quadcol,  	/**< indices for the quadratic part */
		const double *      quadval 	/**< nonzero coefficient of the quadratic part */ );
	
	/** print quadratic constraints */
	DSP_RTN_CODE printQuadRows (const int s);

	/** print quadratic constraints for given params */
	DSP_RTN_CODE printQuadRows (QcRowDataScen *qcdata);

	/** get parameters for CPX interface */
	QcRowDataScen * getQRowsCPXParams (int s) {return &QcRowData_[s];}

	/** get file name */
	const char * getFileName() const {return filename_;} /* must be the same as getNumScenarios() in StoModel */
	
	/** get number of scenarios */
	int getNumScens() const {return nscen_;} /* must be the same as getNumScenarios() in StoModel */
	
	/** get number of quadratic constraints */
	int getNumQRows(int scen) {return QcRowData_[scen].nqrows_;}

	/** get number of non-zero linear terms of quadratic constraints */
	int * getLinearNonZeroCounts(int scen) const {return QcRowData_[scen].linnzcnt_;}

	/** get number of non-zero quadratic terms of quadratic constraints */
	int * getQuadNonZeroCounts(int scen) const {return QcRowData_[scen].quadnzcnt_;}

	/** get rhs of quadratic constraints */
	double * getRhs(int scen) const {return QcRowData_[scen].rhs_;}

	/** get sense of quadratic constraints */
	int * getSense(int scen) const {return QcRowData_[scen].sense_;}

	/** get indices for linear terms of quadratic constraints */
	const int ** getLinearIndices(int scen) const {return (const int **) QcRowData_[scen].linind_;}

	/** get coefficients for linear terms of quadratic constraints */
	const double ** getLinearVals(int scen) const {return (const double **) QcRowData_[scen].linval_;}

	/** get first indices for quadratic terms of quadratic constraints */
	const int ** getQuadIndices1st(int scen) const {return (const int **) QcRowData_[scen].quadrow_;}

	/** get second indices for quadratic terms of quadratic constraints */
	const int ** getQuadIndices2nd(int scen) const {return (const int **) QcRowData_[scen].quadcol_;}

	/** get coefficients for quadratic terms of quadratic constraints */
	const double ** getQuadraticVals(int scen) const {return (const double **) QcRowData_[scen].quadval_;}

protected:

	int				nscen_;
	const char * filename_;

	vector<QcRowDataScen> QcRowData_;
};

#endif /* QCMODEL_H_ */
