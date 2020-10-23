/*
 * QCModel.cpp
 *
 *  Created on: Oct 22, 2020
 *      Author: geunyeong byeon
 */

//#define DSP_DEBUG

#include "Utility/DspMessage.h"
#include "Model/QcModel.h"

QcModel::QcModel(){}

/** copy constructor */
QcModel::QcModel(const QcModel & rhs) {}

QcModel::~QcModel()
{
	FREE_2D_ARRAY_PTR(nscen_, linnzcnt_);
	FREE_2D_ARRAY_PTR(nscen_, quadnzcnt_);
	FREE_2D_ARRAY_PTR(nscen_, rhs_);
	FREE_2D_ARRAY_PTR(nscen_, sense_);

	for (int i = 0; i < nscen_; i++)
	{
		FREE_2D_ARRAY_PTR(nqconstrs_[i], linind_);
		FREE_2D_ARRAY_PTR(nqconstrs_[i], linval_);
		FREE_2D_ARRAY_PTR(nqconstrs_[i], quadrow_);
		FREE_2D_ARRAY_PTR(nqconstrs_[i], quadcol_);
		FREE_2D_ARRAY_PTR(nqconstrs_[i], quadval_);
	}
	FREE_ARRAY_PTR(linind_);
	FREE_ARRAY_PTR(linval_);
	FREE_ARRAY_PTR(quadrow_);
	FREE_ARRAY_PTR(quadcol_);
	FREE_ARRAY_PTR(quadval_);
	FREE_ARRAY_PTR(nqconstrs_);
}

/** set dimensions for second-stage quadratic constraints */
DSP_RTN_CODE QcModel::setQuadDimensions(int nscen, int * nqconstrs)
{
	BGN_TRY_CATCH

	nscen_ = nscen;

	nqconstrs_ = new int [nscen];
	linnzcnt_ = new int * [nscen];
	quadnzcnt_ = new int * [nscen];
	rhs_ = new double * [nscen];
	sense_ = new int * [nscen];

	linind_ = new int ** [nscen];
	linval_ = new double ** [nscen];

	quadrow_ = new int ** [nscen];
	quadcol_ = new int ** [nscen];
	quadval_ = new double ** [nscen];

	for (int s = 0; s < nscen; s++) 
	{
		nqconstrs_[s] = nqconstrs[s]; 	
		linnzcnt_[s] = new int [nqconstrs_[s]];
		quadnzcnt_[s] = new int [nqconstrs_[s]];
		rhs_[s] = new double [nqconstrs_[s]];
		sense_[s] = new int [nqconstrs_[s]];

		linind_[s] = new int * [nqconstrs_[s]];
		linval_[s] = new double * [nqconstrs_[s]];
		quadrow_[s] = new int * [nqconstrs_[s]];
		quadcol_[s] = new int * [nqconstrs_[s]];
		quadval_[s] = new double * [nqconstrs_[s]];
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	
	return DSP_RTN_OK;
}

/** load quadratic constraints to the second stage */
DSP_RTN_CODE QcModel::loadQuadraticConstrs(
        const int           s,     		/**< scenario index */
		const int 			nqconstrs,
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
		const double *      quadval 	/**< nonzero coefficient of the quadratic part */ ){

    BGN_TRY_CATCH

	assert(nqconstrs == nqconstrs_[s]);

	/** allocate values */
	for (int k = 0; k < nqconstrs; k++) 
	{
		linnzcnt_[s][k] = linnzcnt[k];
		quadnzcnt_[s][k] = quadnzcnt[k];
		rhs_[s][k] = rhs[k];
		sense_[s][k] = sense[k];

		linind_[s][k] = new int [linnzcnt[k]];
		linval_[s][k] = new double [linnzcnt[k]];

		quadrow_[s][k] = new int [quadnzcnt[k]];
		quadcol_[s][k] = new int [quadnzcnt[k]];
		quadval_[s][k] = new double [quadnzcnt[k]];

		for (int t = 0; t < linnzcnt_[s][k]; t++) 
		{
			linind_[s][k][t] = linind[linstart[k] + t];
			linval_[s][k][t] = linval[linstart[k] + t];
		}
		
		for (int t = 0; t < quadnzcnt[k]; t++) 
		{
			quadrow_[s][k][t] = quadrow[quadstart[k] + t];
			quadcol_[s][k][t] = quadcol[quadstart[k] + t];
			quadval_[s][k][t] = quadval[quadstart[k] + t];
		}
	}

    END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE QcModel::printQuadraticConstrs (const int s)
{
	BGN_TRY_CATCH

	for (int i = 0; i < nqconstrs_[s]; i++) 
	{
		cout << "Scen " << s << "th " << i << "th quad constr: ";

		for (int lt = 0; lt < linnzcnt_[s][i]; lt++)
		{
			cout << linval_[s][i][lt] << " x" << linind_[s][i][lt] << " + ";
		}
		for (int qt = 0; qt < quadnzcnt_[s][i]-1; qt++)
		{
			cout << quadval_[s][i][qt] << " x" << quadrow_[s][i][qt] << " x" << quadcol_[s][i][qt] << " + ";
		}
		cout << quadval_[s][i][quadnzcnt_[s][i]-1] << " x" << quadrow_[s][i][quadnzcnt_[s][i]-1] << " x" << quadcol_[s][i][quadnzcnt_[s][i]-1];
		if (sense_[s][i] == 'L')
			cout << " <= " << rhs_[s][i] << endl;
		else 
			cout << " >= " << rhs_[s][i] << endl;

	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}


DSP_RTN_CODE QcModel::getQConstrParametersCPX (int s, int &nqconstrs, int *linnzcnt, int * quadnzcnt, double * rhs, int * sense, const int ** linind, const double ** linval, const int ** quadrow, const int ** quadcol, const double ** quadval)
{
	BGN_TRY_CATCH
	
	nqconstrs = getNumQConstrs(s);
	linnzcnt =  getLinearNonZeroCounts(s);
	quadnzcnt = getQuadNonZeroCounts(s); 
	rhs = getRhs(s); 
	sense = getSense(s); 
	linind = getLinearIndices(s); 
	linval = getLinearVals(s); 
	quadrow = getQuadIndices1st(s); 
	quadcol = getQuadIndices2nd(s); 
	quadval = getQuadraticVals(s); 
	
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
