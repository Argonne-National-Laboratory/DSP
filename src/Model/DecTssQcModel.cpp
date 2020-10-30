/*
 * DecTssQcModel.cpp
 *
 *  Created on: Oct 27, 2020
 *      Author: geunyeongbyeon
 */

#include "Utility/DspMessage.h"
#include "Model/DecTssQcModel.h"

DecTssQcModel::DecTssQcModel() :
DecTssModel() {
	/** nothing to do */
}

/** copy constructor */
DecTssQcModel::DecTssQcModel(const DecTssQcModel & rhs) :
DecTssModel(rhs) {
	/** nothing to do */
}

/** copy constructor */
DecTssQcModel::DecTssQcModel(const DecTssModel & rhs) :
DecTssModel(rhs) {
	/** nothing to do */
}
DecTssQcModel::~DecTssQcModel() {
	
	for (int s = 0; s < nscen_; s++)
	{
		FREE_ARRAY_PTR(QcRowData_[s].linnzcnt_);
		FREE_ARRAY_PTR(QcRowData_[s].quadnzcnt_);
		FREE_ARRAY_PTR(QcRowData_[s].rhs_);
		FREE_ARRAY_PTR(QcRowData_[s].sense_);

		FREE_2D_ARRAY_PTR(QcRowData_[s].nqrows_, QcRowData_[s].linind_);
		FREE_2D_ARRAY_PTR(QcRowData_[s].nqrows_, QcRowData_[s].linval_);
		FREE_2D_ARRAY_PTR(QcRowData_[s].nqrows_, QcRowData_[s].quadrow_);
		FREE_2D_ARRAY_PTR(QcRowData_[s].nqrows_, QcRowData_[s].quadcol_);
		FREE_2D_ARRAY_PTR(QcRowData_[s].nqrows_, QcRowData_[s].quadval_);
	}
}


/** construct a map that maps variable names to their indices */
bool DecTssQcModel::mapVarnameIndex(map<string, int> &map_varName_index, const char * smps) 
{
	char core[128];
	sprintf(core, "%s.cor", smps); 
	ifstream corefile (core);

	string name, item;

	map_varName_index.clear();
    int nvars = 0;
	vector<string> rowNames;

	map<string, int>::iterator map_it;
    vector<string>::iterator vec_it;

    bool foundColumns = false;

    if (corefile.is_open()) {
        
	    while (corefile >> item) {
		
            if (item.find("ROWS") != string::npos) {
			
                while (1) {
			        corefile >> item >> name;

                    if (item.find("COLUMNS") != string::npos) {
                        foundColumns = true;
                        map_varName_index[name] = nvars;
                        nvars++;
                        break;
                    }

                    if (name == "ENDATA" || item == "ENDATA") {
                        cout << "Encountered ENDATA before reading Columns" << endl;
                        break;
                    }

			        rowNames.push_back(name);
                }
		
                if (foundColumns) {
                
                    while(1) {
				
	        			corefile >> name >> item;
				
			        	vec_it = find(rowNames.begin(), rowNames.end(), name);
				        if (vec_it == rowNames.end()) 
				        {
					        /** var term */
					        corefile >> item;

					        if (name.find("RHS") != string::npos)
						        break;
					
					        map_it = map_varName_index.find(name);
  					        if (map_it == map_varName_index.end())
    					    {
    					        map_varName_index[name] = nvars;
                                nvars++;
					        }
				        } else 
				        {
					        /* constraint term
					        * does not do anything */
				        }
			        }
		        } else {
                    cout << "Columns data should follow Rows data" << endl;
                    break;
                }
            }

		if (item.find("ENDATA") != string::npos) 
			break;
        
	    }
    corefile.close();
    } else {
        cout << "Unable to open core file";
        return false;
    }

#ifdef DSP_DEBUG
    cout << "variable, index" << endl;
    for (auto &v : map_varName_index)
    cout << v.first << ", " << v.second << endl;
#endif              

	return true;
}
/** read quadratic data file */
DSP_RTN_CODE DecTssQcModel::readQuad(const char * smps, const char * filename)
{
	BGN_TRY_CATCH

	map<string, int> map_varName_index;

	if (!mapVarnameIndex(map_varName_index, smps)) 
	{
		char msg[128];
		sprintf(msg, "Unable to map variables to their indices\n");
		throw msg;
	}

	setFileName(filename);
	char quadfile[128];
	sprintf(quadfile, "%s.txt", filename); 
	ifstream myfile(quadfile);

	string item;
	string name, name2; 
	double val;
	int i, j;
	int sind;
	bool end_of_scen_data = false;

	map<string, int> map_qrowName_index;
	map<string, int>::iterator it;
	vector<char> qrow_sense;

	/* allocate memory */
	setQcDimensions();

	if (myfile.is_open()) {
		while (myfile >> item) {
		
			if (item.find("NAME") != string::npos) 
			{
				myfile >> item;
			} 
			else if (item.find("SCEN") != string::npos) 
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
					QcRowData_[sind].nqrows_ = 0;

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
						map_qrowName_index[name] = QcRowData_[sind].nqrows_;
						QcRowData_[sind].nqrows_++;
					}

					/** allocate memory */
					assert(QcRowData_[sind].nqrows_ == qrow_sense.size());
					assert(QcRowData_[sind].nqrows_ == map_qrowName_index.size());
				
					setQcDimensions(sind, QcRowData_[sind].nqrows_);

				/** read sense_ */
				for (int i = 0; i < QcRowData_[sind].nqrows_; i++) {
					QcRowData_[sind].sense_[i] = qrow_sense[i];
				}

				/** read linind_, linval_ */
				vector<vector<int>> linind(QcRowData_[sind].nqrows_);
				vector<vector<double>> linval(QcRowData_[sind].nqrows_);
				
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
			
				for (i = 0; i < QcRowData_[sind].nqrows_; i++)
				{
					assert(linind[i].size() == linval[i].size());
					QcRowData_[sind].linnzcnt_[i] = linval[i].size();
					QcRowData_[sind].linind_[i] = new int [QcRowData_[sind].linnzcnt_[i]];
					QcRowData_[sind].linval_[i] = new double [QcRowData_[sind].linnzcnt_[i]];

					for (j = 0; j < QcRowData_[sind].linnzcnt_[i]; j++) 
					{
						QcRowData_[sind].linind_[i][j] = linind[i][j];
						QcRowData_[sind].linval_[i][j] = linval[i][j];
					}
				}

				/** read quadrow_, quadcol_, quadval_ */
				vector<vector<int>> quadrow(QcRowData_[sind].nqrows_);
				vector<vector<int>> quadcol(QcRowData_[sind].nqrows_);
				vector<vector<double>> quadval(QcRowData_[sind].nqrows_);	

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

				for (i = 0; i < QcRowData_[sind].nqrows_; i++)
				{
					assert(quadrow[i].size() == quadcol[i].size());
					assert(quadrow[i].size() == quadval[i].size());
					QcRowData_[sind].quadnzcnt_[i] = quadrow[i].size();
					QcRowData_[sind].quadrow_[i] = new int [QcRowData_[sind].quadnzcnt_[i]];
					QcRowData_[sind].quadcol_[i] = new int [QcRowData_[sind].quadnzcnt_[i]];
					QcRowData_[sind].quadval_[i] = new double [QcRowData_[sind].quadnzcnt_[i]];

					for (j = 0; j < QcRowData_[sind].quadnzcnt_[i]; j++) 
					{
						QcRowData_[sind].quadrow_[i][j] = quadrow[i][j];
						QcRowData_[sind].quadcol_[i][j] = quadcol[i][j];
						QcRowData_[sind].quadval_[i][j] = quadval[i][j];
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
							QcRowData_[sind].rhs_[map_qrowName_index[item]] = val;
						}
					}

					if (end_of_scen_data)
						break;
				}
			} else 
			{
				char msg[128];
				sprintf(msg, "Quadratic data file encountered unexpected input: %s\n", item.c_str());
				throw msg;
			}
		}
	    myfile.close();
    } else {
        cout << "Unable to open quad file";
        return 1;
    }

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** set dimensions for second-stage quadratic constraints */
DSP_RTN_CODE DecTssQcModel::setQcDimensions(int * nqrows)
{
	BGN_TRY_CATCH

	QcRowData_.resize(nscen_);

	for (int s = 0; s < nscen_; s++) 
	{
		setQcDimensions(s, nqrows[s]);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	
	return DSP_RTN_OK;
}

DSP_RTN_CODE DecTssQcModel::setQcDimensions(int s, int nqrows)
{
	BGN_TRY_CATCH

	QcRowData_[s].nqrows_ = nqrows; 	
	
	QcRowData_[s].linnzcnt_ = new int [nqrows];
	QcRowData_[s].quadnzcnt_ = new int [nqrows];
	QcRowData_[s].rhs_ = new double [nqrows];
	QcRowData_[s].sense_ = new int [nqrows];

	QcRowData_[s].linind_ = new int * [nqrows];
	QcRowData_[s].linval_ = new double * [nqrows];
	QcRowData_[s].quadrow_ = new int * [nqrows];
	QcRowData_[s].quadcol_ = new int * [nqrows];
	QcRowData_[s].quadval_ = new double * [nqrows];

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	
	return DSP_RTN_OK;
}

/** load quadratic constraints to the second stage */
DSP_RTN_CODE DecTssQcModel::loadQuadraticRows(
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
		const double *      quadval 	/**< nonzero coefficient of the quadratic part */ ){

    BGN_TRY_CATCH

	assert(nqrows == nqrows_[s]);
	
	/** allocate values */
	for (int k = 0; k < nqrows; k++) 
	{
		QcRowData_[s].linnzcnt_[k] = linnzcnt[k];
		QcRowData_[s].quadnzcnt_[k] = quadnzcnt[k];
		QcRowData_[s].rhs_[k] = rhs[k];
		QcRowData_[s].sense_[k] = sense[k];

		QcRowData_[s].linind_[k] = new int [linnzcnt[k]];
		QcRowData_[s].linval_[k] = new double [linnzcnt[k]];

		QcRowData_[s].quadrow_[k] = new int [quadnzcnt[k]];
		QcRowData_[s].quadcol_[k] = new int [quadnzcnt[k]];
		QcRowData_[s].quadval_[k] = new double [quadnzcnt[k]];

		for (int t = 0; t < QcRowData_[s].linnzcnt_[k]; t++) 
		{
			QcRowData_[s].linind_[k][t] = linind[linstart[k] + t];
			QcRowData_[s].linval_[k][t] = linval[linstart[k] + t];
		}
		
		for (int t = 0; t < QcRowData_[s].quadnzcnt_[k]; t++) 
		{
			QcRowData_[s].quadrow_[k][t] = quadrow[quadstart[k] + t];
			QcRowData_[s].quadcol_[k][t] = quadcol[quadstart[k] + t];
			QcRowData_[s].quadval_[k][t] = quadval[quadstart[k] + t];
		}
	}

    END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DecTssQcModel::printQuadRows (const int s)
{
	BGN_TRY_CATCH

	for (int i = 0; i < QcRowData_[s].nqrows_; i++) 
	{
		cout << "Scen " << s << "th " << i << "th quad constr: ";

		for (int lt = 0; lt < QcRowData_[s].linnzcnt_[i]; lt++)
		{
			cout << QcRowData_[s].linval_[i][lt] << " x" << QcRowData_[s].linind_[i][lt] << " + ";
		}
		for (int qt = 0; qt < QcRowData_[s].quadnzcnt_[i]-1; qt++)
		{
			cout << QcRowData_[s].quadval_[i][qt] << " x" << QcRowData_[s].quadrow_[i][qt] << " x" << QcRowData_[s].quadcol_[i][qt] << " + ";
		}
		cout << QcRowData_[s].quadval_[i][QcRowData_[s].quadnzcnt_[i]-1] << " x" << QcRowData_[s].quadrow_[i][QcRowData_[s].quadnzcnt_[i]-1] << " x" << QcRowData_[s].quadcol_[i][QcRowData_[s].quadnzcnt_[i]-1];
		if (QcRowData_[s].sense_[i] == 'L')
			cout << " <= " << QcRowData_[s].rhs_[i] << endl;
		else 
			cout << " >= " << QcRowData_[s].rhs_[i] << endl;

	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DecTssQcModel::printQuadRows (QcRowDataScen *qcdata)
{
	BGN_TRY_CATCH

	for (int i = 0; i < qcdata->nqrows_; i++) 
	{
		cout << i << "th quad constr: ";

		for (int lt = 0; lt < qcdata->linnzcnt_[i]; lt++)
		{
			cout << qcdata->linval_[i][lt] << " x" << qcdata->linind_[i][lt] << " + ";
		}
		for (int qt = 0; qt < qcdata->quadnzcnt_[i]-1; qt++)
		{
			cout << qcdata->quadval_[i][qt] << " x" << qcdata->quadrow_[i][qt] << " x" << qcdata->quadcol_[i][qt] << " + ";
		}
		cout << qcdata->quadval_[i][qcdata->quadnzcnt_[i]-1] << " x" << qcdata->quadrow_[i][qcdata->quadnzcnt_[i]-1] << " x" << qcdata->quadcol_[i][qcdata->quadnzcnt_[i]-1];
		if (qcdata->sense_[i] == 'L')
			cout << " <= " << qcdata->rhs_[i] << endl;
		else 
			cout << " >= " << qcdata->rhs_[i] << endl;

	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
