/*
 * QCModel.cpp
 *
 *  Created on: Oct 22, 2020
 *      Author: geunyeong byeon
 */

//#define DSP_DEBUG

#include "Utility/DspMessage.h"
#include "Model/QcModel.h"

QcModel::QcModel(){
	nscen_ = -1;
}

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

// /** construct a map that maps variable names to their indices */
// DSP_RTN_CODE mapVarnameIndex(vector<string> &varNames, const char * corefilename) 
// {
// 	BGN_TRY_CATCH

// 	ifstream corefile (corefilename, ios::in | ios:: binary);

// 	string name, item, rowname, varname;

// 	varNames.clear();
// 	vector<string> rowNames;

// 	std::vector<string>::iterator it;

// 	while (corefile >> item) {
		
// 		if (item.find("ROWS") != string::npos) {
			
// 			corefile >> name >> item;
// 			rowNames.push_back(name);
// 		}

// 		if (item.find("COLUMNS") != string::npos) {
			
// 			corefile >> varname;
// 			varNames.push_back(varname);

// 			while (1) {
				
// 				corefile >> name >> item;
				
// 				it = find(rowNames.begin(), rowNames.end(), name);
// 				if (it == rowNames.end()) 
// 				{
// 					/** var term */
// 					corefile >> item;

// 					if (name == "RHS")
// 						break;
					
// 					it = find(varNames.begin(), varNames.end(), name);
//   					if (it != varNames.end())
//     					std::cout << "Element found in myvector: " << *it << '\n';
//   					else {
//     					std::cout << "Element not found in myvector\n";
// 						varNames.push_back(name);
// 					}
// 				} else 
// 				{
// 					/* constraint term
// 					 * does not do anything */
// 				}
// 			}
// 		}

// 		if (item.find("ENDATA") != string::npos) 
// 			break;
// 	}

// 	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

// 	return DSP_RTN_OK;

// }
/** read quadratic data file */
DSP_RTN_CODE QcModel::readQuad(map<string, int> &map_varName_index, const char * filename) 
{
	BGN_TRY_CATCH

	ifstream myfile(filename);

	string item;
	string name, name2; 
	double val;
	int i, j;
	int sind;
	bool end_of_scen_data = false;

	map<string, int> map_qrowName_index;
	map<string, int>::iterator it;
	vector<char> qrow_sense;

	while (myfile >> item) {
		
		if (item.find("NAME") != string::npos) 
		{
			myfile >> item;
		} else if (item.find("NSCEN") != string::npos) 
		{
			int nscen;
			myfile >> nscen;

			/** allocate memory */
			setQuadDimensions(nscen);
		} else if (item.find("SCEN") != string::npos) 
		{	
			/** start reading quad constr data for some scenario */		
			if (nscen_ < 0)
			{
				char msg[128];
				sprintf(msg, "NSCEN should be provided before SCEN\n");
				throw msg;
			}

			while (1) {
			myfile >> sind;
			
			nqconstrs_[sind] = 0;

			myfile >> item;
			if (item.find("QUADROWS") == string::npos) {
				char msg[128];
				sprintf(msg, "Quadratic Constraints Data for each SCEN should be provided in this order: QUADROWS, LINTERMS, QUADTERMS, RHS\n");
				throw msg;
			} 
			
			/** read row data */
			while (1) 
			{
				myfile >> item;

				if (item.find("LINTERMS") != string::npos)
					break;
				else if (item == "G")
					qrow_sense.push_back('G');
				else if (item == "L")
					qrow_sense.push_back('L');
				else {
					char msg[128];
					sprintf(msg, "Quadratic constraints sense must be 'G' or 'L'\n");
					throw msg;
				} 

				myfile >> name;

				map_qrowName_index[name] = nqconstrs_[sind];

				nqconstrs_[sind]++;
			}

			/** allocate memory */
			assert(nqconstrs_[sind] == qrow_sense.size());
			assert(nqconstrs_[sind] == map_qrowName_index.size());
			
			setQuadDimensions(sind, nqconstrs_[sind]);

			/** read sense_ */
			for (int i = 0; i < nqconstrs_[sind]; i++) {
				sense_[sind][i] = qrow_sense[i];
			}

			/** read linind_, linval_ */
			vector<vector<int>> linind(nqconstrs_[sind]);
			vector<vector<double>> linval(nqconstrs_[sind]);
				
			while (1) {
				myfile >> item;

				if (item.find("QUADTERMS") != string::npos)
					break;
				
				it = map_qrowName_index.find(item);

				if (it == map_qrowName_index.end()) 
				{
					char msg[128];
					sprintf(msg, "All rows should be declared before SCEN data\n");
					throw msg; 
				} else 
				{
					myfile >> name >> val;
					linind[map_qrowName_index[item]].push_back(map_varName_index[name]);
					linval[map_qrowName_index[item]].push_back(val);
				}
			}
			
			for (i = 0; i < nqconstrs_[sind]; i++)
			{
				assert(linind[i].size() == linval[i].size());
				linnzcnt_[sind][i] = linval[i].size();
				linind_[sind][i] = new int [linnzcnt_[sind][i]];
				linval_[sind][i] = new double [linnzcnt_[sind][i]];

				for (j = 0; j < linnzcnt_[sind][i]; j++) 
				{
					linind_[sind][i][j] = linind[i][j];
					linval_[sind][i][j] = linval[i][j];
				}
			}

			/** read quadrow_, quadcol_, quadval_ */
			vector<vector<int>> quadrow(nqconstrs_[sind]);
			vector<vector<int>> quadcol(nqconstrs_[sind]);
			vector<vector<double>> quadval(nqconstrs_[sind]);

			while (1) {
				myfile >> item;

				if (item.find("RHS") != string::npos)
					break;

				it = map_qrowName_index.find(item);

				if (it == map_qrowName_index.end()) 
				{
					char msg[128];
					sprintf(msg, "All rows should be declared before SCEN data\n");
					throw msg; 
				} else 
				{
					myfile >> name >> name2 >> val;
					quadrow[map_qrowName_index[item]].push_back(map_varName_index[name]);
					quadcol[map_qrowName_index[item]].push_back(map_varName_index[name2]);
					quadval[map_qrowName_index[item]].push_back(val);
				}

			}

			for (i = 0; i < nqconstrs_[sind]; i++)
			{
				assert(quadrow[i].size() == quadcol[i].size());
				assert(quadrow[i].size() == quadval[i].size());
				quadnzcnt_[sind][i] = quadrow[i].size();
				quadrow_[sind][i] = new int [quadnzcnt_[sind][i]];
				quadcol_[sind][i] = new int [quadnzcnt_[sind][i]];
				quadval_[sind][i] = new double [quadnzcnt_[sind][i]];

				for (j = 0; j < quadnzcnt_[sind][i]; j++) 
				{
					quadrow_[sind][i][j] = quadrow[i][j];
					quadcol_[sind][i][j] = quadcol[i][j];
					quadval_[sind][i][j] = quadval[i][j];
				}
			}

			/** read rhs_ */
			while (1) {
				myfile >> item;

				if (item.find("SCEN") != string::npos)
					break;
				else if (item.find("ENDATA") != string::npos) 
				{
					end_of_scen_data = true;
					break;
				}
			
				it = map_qrowName_index.find(item);

					if (it == map_qrowName_index.end()) 
					{
						char msg[128];
						sprintf(msg, "All rows should be declared before SCEN data\n");
						throw msg; 
					} else 
					{
						myfile >> val;
						rhs_[sind][map_qrowName_index[item]] = val;
					}
				}

				if (end_of_scen_data)
					break;

			}
		} 
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** set dimensions for second-stage quadratic constraints */
DSP_RTN_CODE QcModel::setQuadDimensions(int nscen, int * nqconstrs)
{
	BGN_TRY_CATCH

	setQuadDimensions(nscen);

	for (int s = 0; s < nscen; s++) 
	{
		setQuadDimensions(s, nqconstrs[s]);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	
	return DSP_RTN_OK;
}
DSP_RTN_CODE QcModel::setQuadDimensions(int nscen)
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

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	
	return DSP_RTN_OK;
}
DSP_RTN_CODE QcModel::setQuadDimensions(int s, int nqconstrs)
{
	BGN_TRY_CATCH

	nqconstrs_[s] = nqconstrs; 	
	
	linnzcnt_[s] = new int [nqconstrs];
	quadnzcnt_[s] = new int [nqconstrs];
	rhs_[s] = new double [nqconstrs];
	sense_[s] = new int [nqconstrs];

	linind_[s] = new int * [nqconstrs];
	linval_[s] = new double * [nqconstrs];
	quadrow_[s] = new int * [nqconstrs];
	quadcol_[s] = new int * [nqconstrs];
	quadval_[s] = new double * [nqconstrs];

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
