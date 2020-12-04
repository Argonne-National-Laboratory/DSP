/*
 * StoModel.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: kibaekkim
 */
// #define DSP_DEBUG
#include "CoinHelperFunctions.hpp"
#include "StoModel.h"

StoModel::StoModel() :
		nscen_(0),
		nstgs_(0),
		nrows_(NULL),
		ncols_(NULL),
		nints_(NULL),
		rstart_(NULL),
		cstart_(NULL),
		nrows_core_(0),
		ncols_core_(0),
		nints_core_(0),
		rows_core_(NULL),
		clbd_core_(NULL),
		cubd_core_(NULL),
		obj_core_(NULL),
		qobj_core_(NULL),
		rlbd_core_(NULL),
		rubd_core_(NULL),
		ctype_core_(NULL),
		prob_(NULL),
		mat_scen_(NULL),
		clbd_scen_(NULL),
		cubd_scen_(NULL),
		obj_scen_(NULL),
		qobj_scen_(NULL),
		rlbd_scen_(NULL),
		rubd_scen_(NULL),
		qc_row_core_(NULL),
		qc_row_scen_(NULL),
		fromSMPS_(false),
		isdro_(false),
		nrefs_(0),
		wass_eps_(0.0),
		wass_dist_(NULL),
		refs_probability_(NULL)
{
	/** nothing to do */
}

/** copy constructor */
StoModel::StoModel(const StoModel & rhs) :
		nscen_(rhs.nscen_),
		nstgs_(rhs.nstgs_),
		nrows_core_(rhs.nrows_core_),
		ncols_core_(rhs.ncols_core_),
		nints_core_(rhs.nints_core_),
		fromSMPS_(rhs.fromSMPS_),
		isdro_(rhs.isdro_),
		nrefs_(rhs.nrefs_),
		wass_eps_(rhs.wass_eps_),
		wass_dist_(rhs.wass_dist_),
		refs_probability_(rhs.refs_probability_),
		qc_row_core_(rhs.qc_row_core_)
{
	/** allocate memory */
	nrows_      = new int [nstgs_];
	ncols_      = new int [nstgs_];
	nints_      = new int [nstgs_];
	rstart_     = new int [nstgs_];
	cstart_     = new int [nstgs_];
	rows_core_  = new CoinPackedVector * [nrows_core_];
	clbd_core_  = new double * [nstgs_];
	cubd_core_  = new double * [nstgs_];
	obj_core_   = new double * [nstgs_];
	qobj_core_  = new CoinPackedMatrix * [nstgs_];
	rlbd_core_  = new double * [nstgs_];
	rubd_core_  = new double * [nstgs_];
	ctype_core_ = new char * [nstgs_];
	prob_       = new double [nscen_];
	mat_scen_   = new CoinPackedMatrix * [nscen_];
	clbd_scen_  = new CoinPackedVector * [nscen_];
	cubd_scen_  = new CoinPackedVector * [nscen_];
	obj_scen_   = new CoinPackedVector * [nscen_];
	qobj_scen_  = new CoinPackedMatrix * [nscen_];
	rlbd_scen_  = new CoinPackedVector * [nscen_];
	rubd_scen_  = new CoinPackedVector * [nscen_];
	qc_row_scen_ = new QuadRowData [nscen_];

	/** copy */
	CoinCopyN(rhs.nrows_, nstgs_, nrows_);
	CoinCopyN(rhs.ncols_, nstgs_, ncols_);
	CoinCopyN(rhs.nints_, nstgs_, nints_);
	CoinCopyN(rhs.rstart_, nstgs_, rstart_);
	CoinCopyN(rhs.cstart_, nstgs_, cstart_);
	for (int i = nrows_core_ - 1; i >= 0; --i)
	{
		if (rhs.rows_core_[i])
			rows_core_[i] = new CoinPackedVector(*(rhs.rows_core_[i]));
		else
			rows_core_[i] = new CoinPackedVector;
	}
	for (int i = nstgs_ - 1; i >= 0; --i)
	{
		clbd_core_[i] = new double [ncols_[i]];
		cubd_core_[i] = new double [ncols_[i]];
		obj_core_[i] = new double [ncols_[i]];
		/** copy quadratic objective information */
		if (rhs.qobj_core_[i]!=NULL){
			qobj_core_[i] = new CoinPackedMatrix(*(rhs.qobj_core_[i]));
		}
		ctype_core_[i] = new char [ncols_[i]];
		rlbd_core_[i] = new double [nrows_[i]];
		rubd_core_[i] = new double [nrows_[i]];
		CoinCopyN(rhs.clbd_core_[i], ncols_[i], clbd_core_[i]);
		CoinCopyN(rhs.cubd_core_[i], ncols_[i], cubd_core_[i]);
		CoinCopyN(rhs.obj_core_[i], ncols_[i], obj_core_[i]);
		CoinCopyN(rhs.rlbd_core_[i], nrows_[i], rlbd_core_[i]);
		CoinCopyN(rhs.rubd_core_[i], nrows_[i], rubd_core_[i]);
		CoinCopyN(rhs.ctype_core_[i], ncols_[i], ctype_core_[i]);
	}
	for (int i = nscen_ - 1; i >= 0; --i)
	{
		qc_row_scen_[i] = rhs.qc_row_scen_[i];
		prob_[i] = rhs.prob_[i];
		if (rhs.mat_scen_[i])
			mat_scen_[i] = new CoinPackedMatrix(*(rhs.mat_scen_[i]));
		else
			mat_scen_[i] = new CoinPackedMatrix;
		if (rhs.clbd_scen_[i])
			clbd_scen_[i] = new CoinPackedVector(*(rhs.clbd_scen_[i]));
		else
			clbd_scen_[i] = new CoinPackedVector;
		if (rhs.cubd_scen_[i])
			cubd_scen_[i] = new CoinPackedVector(*(rhs.cubd_scen_[i]));
		else
			cubd_scen_[i] = new CoinPackedVector;
		if (rhs.obj_scen_[i])
			obj_scen_[i] = new CoinPackedVector(*(rhs.obj_scen_[i]));
		else
			obj_scen_[i] = new CoinPackedVector;
		if (rhs.qobj_scen_[i])
			qobj_scen_[i] = new CoinPackedMatrix(*(rhs.qobj_scen_[i]));
		else
			qobj_scen_[i] = new CoinPackedMatrix;
		if (rhs.rlbd_scen_[i])
			rlbd_scen_[i] = new CoinPackedVector(*(rhs.rlbd_scen_[i]));
		else
			rlbd_scen_[i] = new CoinPackedVector;
		if (rhs.rubd_scen_[i])
			rubd_scen_[i] = new CoinPackedVector(*(rhs.rubd_scen_[i]));
		else
			rubd_scen_[i] = new CoinPackedVector;
	}

	/** copy initial solutions */
	for (unsigned i = 0; i < rhs.init_solutions_.size(); ++i)
		init_solutions_.push_back(new CoinPackedVector(rhs.init_solutions_[i]));
}

StoModel::~StoModel()
{
	FREE_ARRAY_PTR(nrows_);
	FREE_ARRAY_PTR(ncols_);
	FREE_ARRAY_PTR(nints_);
	FREE_ARRAY_PTR(rstart_);
	FREE_ARRAY_PTR(cstart_);
	FREE_ARRAY_PTR(prob_);
	FREE_ARRAY_PTR(refs_probability_);
	FREE_2D_PTR(nrows_core_, rows_core_);
	FREE_2D_ARRAY_PTR(nstgs_, clbd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, cubd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, obj_core_);
	FREE_2D_PTR(nstgs_, qobj_core_);
	FREE_2D_ARRAY_PTR(nstgs_, rlbd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, rubd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, ctype_core_);
	FREE_2D_ARRAY_PTR(nrefs_, wass_dist_);
	FREE_2D_PTR(nscen_, mat_scen_);
	FREE_2D_PTR(nscen_, clbd_scen_);
	FREE_2D_PTR(nscen_, cubd_scen_);
	FREE_2D_PTR(nscen_, obj_scen_);
	FREE_2D_PTR(nscen_, qobj_scen_);
	FREE_2D_PTR(nscen_, rlbd_scen_);
	FREE_2D_PTR(nscen_, rubd_scen_);
	FREE_PTR(qc_row_core_);
	FREE_ARRAY_PTR(qc_row_scen_);
	scen2stg_.clear();
	nscen_ = 0;
	nstgs_ = 0;
	nrows_core_ = 0;
	ncols_core_ = 0;
	for (unsigned i = 0; i < init_solutions_.size(); ++i)
		FREE_PTR(init_solutions_[i]);
}

DSP_RTN_CODE StoModel::readSmps(const char * filename)
{
	int i, j, s;
	SmiScnModel smi;
	SmiCoreData * core = NULL;
	SmiNodeData * node = NULL;

	BGN_TRY_CATCH

	/** read SMPS files */
	double stime = CoinCpuTime();
	int smistatus = smi.readSmps(filename);
	if (smistatus < 0)
	{
		throw "Failed to read file.";
	}
	std::cout << "Read SMPS files: " << CoinCpuTime() - stime << " sec." << std::endl;
	
	core = smi.getCore();

	/** number of stages */
	nstgs_ = core->getNumStages();
	
	/** Do not support multi-stage */
	if (nstgs_ > 2)
	{
		char tmpstr[128];
		sprintf(tmpstr, "The problem has %d stages. Multi-stage model is not supported.\n", nstgs_);
		CoinError(tmpstr, "readSmps", "StoModel");
	}

	/** number of scenarios */
	nscen_ = smi.getNumScenarios();

	/** allocate memory */
	assert(nstgs_ == 2);
	assert(nscen_ > 0);
	nrows_      = new int [nstgs_];
	ncols_      = new int [nstgs_];
	nints_      = new int [nstgs_];
	rstart_     = new int [nstgs_];
	cstart_     = new int [nstgs_];
	clbd_core_  = new double * [nstgs_];
	cubd_core_  = new double * [nstgs_];
	obj_core_   = new double * [nstgs_];
	qobj_core_  = new CoinPackedMatrix * [nstgs_];
	rlbd_core_  = new double * [nstgs_];
	rubd_core_  = new double * [nstgs_];
	ctype_core_ = new char * [nstgs_];
	prob_       = new double [nscen_];
	mat_scen_   = new CoinPackedMatrix * [nscen_];
	clbd_scen_  = new CoinPackedVector * [nscen_];
	cubd_scen_  = new CoinPackedVector * [nscen_];
	obj_scen_   = new CoinPackedVector * [nscen_];
	qobj_scen_  = new CoinPackedMatrix * [nscen_];
	rlbd_scen_  = new CoinPackedVector * [nscen_];
	rubd_scen_  = new CoinPackedVector * [nscen_];
	qc_row_core_ = new QuadRowData;
	qc_row_scen_ = new QuadRowData [nscen_];

	for (i=0; i<nstgs_;i++){
		qobj_core_[i]=NULL;
	}

	for (i=0; i<nscen_; i++){
		qobj_scen_[i]=NULL;
	}
	/** stage information */
	nrows_core_ = 0;
	ncols_core_ = 0;
	nints_core_ = 0;
	for (i = 0; i < nstgs_; ++i)
	{
		assert(core->getNumRows(i) >= 0);
		assert(core->getNumCols(i) >= 0);
		assert(core->getRowStart(i) >= 0);
		assert(core->getColStart(i) >= 0);
		nrows_[i] = core->getNumRows(i);
		ncols_[i] = core->getNumCols(i);
		nints_[i] = (int) core->getIntCols(i).size();
		rstart_[i] = core->getRowStart(i);
		cstart_[i] = core->getColStart(i);
		nrows_core_ += nrows_[i];
		ncols_core_ += ncols_[i];
		clbd_core_[i] = new double [ncols_[i]];
		cubd_core_[i] = new double [ncols_[i]];
		obj_core_[i] = new double [ncols_[i]];
		rlbd_core_[i] = new double [nrows_[i]];
		rubd_core_[i] = new double [nrows_[i]];
		ctype_core_[i] = new char [ncols_[i]];
		core->copyColLower(clbd_core_[i], i);
		core->copyColUpper(cubd_core_[i], i);
		core->copyObjective(obj_core_[i], i);
		core->copyRowLower(rlbd_core_[i], i);
		core->copyRowUpper(rubd_core_[i], i);

		/** set column types */
		CoinFillN(ctype_core_[i], ncols_[i], 'C');
		for (j = 0; j < core->getBinaryLength(); ++j)
		{
			if (core->getBinaryIndices()[j] < cstart_[i] ||
				core->getBinaryIndices()[j] >= cstart_[i] + ncols_[i])
				continue;
			ctype_core_[i][core->getBinaryIndices()[j] - cstart_[i]] = 'B';
			clbd_core_[i][core->getBinaryIndices()[j] - cstart_[i]] = 0.0;
			cubd_core_[i][core->getBinaryIndices()[j] - cstart_[i]] = 1.0;
			nints_core_++;
		}
		for (j = 0; j < core->getIntegerLength(); ++j)
		{
			if (core->getIntegerIndices()[j] < cstart_[i] ||
				core->getIntegerIndices()[j] >= cstart_[i] + ncols_[i])
				continue;
			ctype_core_[i][core->getIntegerIndices()[j] - cstart_[i]] = 'I';
			nints_core_++;
		}
	}

	/** construct core matrix rows */
	j = 0;
	rows_core_  = new CoinPackedVector * [nrows_core_];
	for (s = 0; s < nstgs_; ++s)
	{
		node = core->getNode(s);
		for (i = rstart_[s]; i < rstart_[s] + nrows_[s]; ++i)
		{
			/**
			 * TODO: I assume the core matrix rows are well ordered.
			 */
			if (i != j)
			{
				CoinError("Unexpected core file structure.", "readSmps", "StoModel");
			}
			rows_core_[j++] = new CoinPackedVector(
					node->getRowLength(i),
					node->getRowIndices(i),
					node->getRowElements(i));
		}
	}
	assert(j == nrows_core_);

	/** clear scenario-stage map */
	scen2stg_.clear();

	std::vector<int> lens;
	for (i = 0; i < nscen_; ++i)
	{
		/** get stage corresponding to scenario */
		int stg = smi.getLeafNode(i)->getStage();

		/** add mapping */
		scen2stg_.insert(std::pair<int,int>(i,stg));

		/** probability */
		prob_[i] = smi.getLeafNode(i)->getProb();

		/** node data */
		node = smi.getLeafNode(i)->getNode();
		clbd_scen_[i] = new CoinPackedVector(
				node->getColLowerLength(), node->getColLowerIndices(), node->getColLowerElements());
		cubd_scen_[i] = new CoinPackedVector(
				node->getColUpperLength(), node->getColUpperIndices(), node->getColUpperElements());
		obj_scen_[i] = new CoinPackedVector(
				node->getObjectiveLength(), node->getObjectiveIndices(), node->getObjectiveElements());
		rlbd_scen_[i] = new CoinPackedVector(
				node->getRowLowerLength(), node->getRowLowerIndices(), node->getRowLowerElements());
		rubd_scen_[i] = new CoinPackedVector(
				node->getRowUpperLength(), node->getRowUpperIndices(), node->getRowUpperElements());
		lens.resize(nrows_[stg]);
		for (j = rstart_[stg]; j < rstart_[stg] + nrows_[stg]; ++j)
			lens[j-rstart_[stg]] = node->getRowLength(j);
		mat_scen_[i]  = new CoinPackedMatrix(false, ncols_[stg], nrows_[stg], node->getNumMatrixElements(), 
				node->getRowElements(rstart_[stg]), node->getRowIndices(rstart_[stg]), node->getRowStarts(rstart_[stg]), &lens[0]);
	}

	/** mark */
	fromSMPS_ = true;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE StoModel::readDro(const char * filename)
{
	if (fromSMPS_ == false) {
		std::cerr << "SMPS should be read first." << std::endl;
		return DSP_RTN_ERR;
	}

	BGN_TRY_CATCH

	string line;
	ifstream myfile(filename);

	if (myfile.is_open()) {
		/**
		 * File contents:
		 * 1. wass_size_
		 * 2. nrefs_
		 * 3. refs_probability_ (size of nrefs by ,)
		 * 4. wass_dist_ (size of nrefs by nscen)
		 */
		int contents_num = 1;
		size_t found;
		string elem;
		while(getline(myfile, line)) {
			// std::cout << line << "," << contents_num << std::endl;
			if (contents_num == 1) {
				wass_eps_ = atof(line.c_str());
				printf("[DRO] The Wasserstein ball size is %e.\n", wass_eps_);
				contents_num++;
			} else if (contents_num == 2) {
				nrefs_ = atoi(line.c_str());
				// printf("nrefs_ has %d.\n", nrefs_);
				printf("[DRO] The first %d scenarios are treated as the reference scenarios.\n", nrefs_);
				contents_num++;
			} else if (contents_num == 3) {
				refs_probability_ = new double [nrefs_];
				
				int i;
				for (i = 0; i < nrefs_; ++i) {
					found = line.find_first_of(" ,");
					// cout << found << endl;
					if(string::npos == found && i < nrefs_ - 1)
						break;
					elem = line.substr(0, found);
					// cout << elem << endl;
					refs_probability_[i] = atof(elem.c_str());
					// printf("refs_probability_[%d] has %f.\n", i, refs_probability_[i]);
					if (string::npos != found)
						line = line.substr(found+1);
					// cout << line << endl;
				}
				if (i != nrefs_) {
					throw "invalid input for reference probability\n";
					break;
				}
				contents_num++;
			} else if (contents_num >= 4) {
				if (contents_num == 4) {
					wass_dist_ = new double * [nrefs_];
					for (int i = 0; i < nrefs_; ++i) {
						wass_dist_[i] = new double [nscen_];
					}
				}

				int i;
				for (i = 0; i < nrefs_; ++i) {
					found = line.find_first_of(" ,");
					if(string::npos == found && i < nrefs_ - 1)
						break;
					wass_dist_[i][contents_num-4] = atof(line.substr(0, found).c_str());
					if (string::npos != found)
						line = line.substr(found+1);
				}
				if (i != nrefs_) {
					throw "invalid number of columns for Wasserstein distance\n";
					break;
				}
				contents_num++;
			}
		}
		if (contents_num-4 != nscen_) {
			throw "invalid number of rows for Wasserstein distance\n";
		}
		isdro_ = true;
	} else {
		char msg[128];
		sprintf(msg, "unable to open file %s\n", filename);
		throw msg;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** read quadratic data file */
DSP_RTN_CODE StoModel::readQuad(const char * smps, const char * filename)
{
	BGN_TRY_CATCH

	map<string, int> map_varName_index;

	if (!mapVarnameIndex(map_varName_index, smps)) 
	{
		char msg[128];
		sprintf(msg, "Unable to map variables to their indices\n");
		throw msg;
	}

	char quadfile[128];
	sprintf(quadfile, "%s.txt", filename); 
	ifstream myfile(quadfile);

	string item;
	string name, name2; 
	double val;
	int i, j;
	int sind;
	bool end_of_scen_data = false;

	map<string, int>::iterator it;

	if (myfile.is_open()) {
		
		while (myfile >> item) {
		
			if (item.find("NAME") != string::npos) 
			{
				myfile >> item;
			} 
			else if (item.find("SCEN") != string::npos) 
			{	
				/** start reading quad constr data for some scenario */		
				while (1) {
					
					myfile >> sind;
					vector<char> qrow_sense;
					map<string, int> map_qrowName_index;

					myfile >> item;
					if (item.find("QUADROWS") == string::npos) {
						char msg[128];
						sprintf(msg, "Quadratic Constraints Data for each SCEN should be provided in this order: QUADROWS, LINTERMS, QUADTERMS, RHS\n");
						throw msg;
					}

					/** read row data */
					int nqrows = 0;
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
						map_qrowName_index[name] = nqrows;
						nqrows++;
					}

					/** allocate memory */
					assert(nqrows == qrow_sense.size());
					assert(nqrows == map_qrowName_index.size());

					/** read linind_, linval_ */
					vector<vector<int>> linind(nqrows);
					vector<vector<double>> linval(nqrows);
					
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
					
					/** read quadrow_, quadcol_, quadval_ */
					vector<vector<int>> quadrow(nqrows);
					vector<vector<int>> quadcol(nqrows);
					vector<vector<double>> quadval(nqrows);	

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

					if (sind < 0) {
						/* if sind == -1, it is qc data for core */
						qc_row_core_->nqrows = nqrows;
						qc_row_core_->linnzcnt = new int [nqrows];
						qc_row_core_->quadnzcnt = new int [nqrows];
						qc_row_core_->rhs = new double [nqrows];
						qc_row_core_->sense = new int [nqrows];
						qc_row_core_->linind = new int * [nqrows];
						qc_row_core_->linval = new double * [nqrows];
						qc_row_core_->quadrow = new int * [nqrows];
						qc_row_core_->quadcol = new int * [nqrows];
						qc_row_core_->quadval = new double * [nqrows];
						
						/** read sense_ */
						for (i = 0; i < nqrows; i++) {
							qc_row_core_->sense[i] = qrow_sense[i];
						
							assert(linind[i].size() == linval[i].size());
							qc_row_core_->linnzcnt[i] = linval[i].size();
							qc_row_core_->linind[i] = new int [qc_row_core_->linnzcnt[i]];
							qc_row_core_->linval[i] = new double [qc_row_core_->linnzcnt[i]];

							assert(quadrow[i].size() == quadcol[i].size());
							assert(quadrow[i].size() == quadval[i].size());
							qc_row_core_->quadnzcnt[i] = quadrow[i].size();
							qc_row_core_->quadrow[i] = new int [qc_row_core_->quadnzcnt[i]];
							qc_row_core_->quadcol[i] = new int [qc_row_core_->quadnzcnt[i]];
							qc_row_core_->quadval[i] = new double [qc_row_core_->quadnzcnt[i]];

							for (j = 0; j < qc_row_core_->linnzcnt[i]; j++) 
							{
								qc_row_core_->linind[i][j] = linind[i][j];
								qc_row_core_->linval[i][j] = linval[i][j];
							}

							for (j = 0; j < qc_row_core_->quadnzcnt[i]; j++) 
							{
								qc_row_core_->quadrow[i][j] = quadrow[i][j];
								qc_row_core_->quadcol[i][j] = quadcol[i][j];
								qc_row_core_->quadval[i][j] = quadval[i][j];
							}
						}
					} else {
						qc_row_scen_[sind].nqrows = nqrows;
						qc_row_scen_[sind].linnzcnt = new int [nqrows];
						qc_row_scen_[sind].quadnzcnt = new int [nqrows];
						qc_row_scen_[sind].rhs = new double [nqrows];
						qc_row_scen_[sind].sense = new int [nqrows];
						qc_row_scen_[sind].linind = new int * [nqrows];
						qc_row_scen_[sind].linval = new double * [nqrows];
						qc_row_scen_[sind].quadrow = new int * [nqrows];
						qc_row_scen_[sind].quadcol = new int * [nqrows];
						qc_row_scen_[sind].quadval = new double * [nqrows];
						
						/** read sense_ */
						for (i = 0; i < nqrows; i++) {
							qc_row_scen_[sind].sense[i] = qrow_sense[i];

							assert(linind[i].size() == linval[i].size());
							qc_row_scen_[sind].linnzcnt[i] = linval[i].size();
							qc_row_scen_[sind].linind[i] = new int [qc_row_scen_[sind].linnzcnt[i]];
							qc_row_scen_[sind].linval[i] = new double [qc_row_scen_[sind].linnzcnt[i]];

							assert(quadrow[i].size() == quadcol[i].size());
							assert(quadrow[i].size() == quadval[i].size());
							qc_row_scen_[sind].quadnzcnt[i] = quadrow[i].size();
							qc_row_scen_[sind].quadrow[i] = new int [qc_row_scen_[sind].quadnzcnt[i]];
							qc_row_scen_[sind].quadcol[i] = new int [qc_row_scen_[sind].quadnzcnt[i]];
							qc_row_scen_[sind].quadval[i] = new double [qc_row_scen_[sind].quadnzcnt[i]];

							for (j = 0; j < qc_row_scen_[sind].linnzcnt[i]; j++) 
							{
								qc_row_scen_[sind].linind[i][j] = linind[i][j];
								qc_row_scen_[sind].linval[i][j] = linval[i][j];
							}

							for (j = 0; j < qc_row_scen_[sind].quadnzcnt[i]; j++) 
							{
								qc_row_scen_[sind].quadrow[i][j] = quadrow[i][j];
								qc_row_scen_[sind].quadcol[i][j] = quadcol[i][j];
								qc_row_scen_[sind].quadval[i][j] = quadval[i][j];
							}
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
								if (sind < 0)
									qc_row_core_->rhs[map_qrowName_index[item]] = val;
								else 
									qc_row_scen_[sind].rhs[map_qrowName_index[item]] = val;
							}
					}
					#ifdef DSP_DEBUG
						printQuadRows(sind);
					#endif
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

/** construct a map that maps variable names to their indices */
bool StoModel::mapVarnameIndex(map<string, int> &map_varName_index, const char * smps) 
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

void StoModel::setSolution(int size, double * solution)
{
	if (size > 0)
		init_solutions_.push_back(new CoinPackedVector(size, solution));
}

DSP_RTN_CODE StoModel::setWassersteinAmbiguitySet(double lp_norm, double eps)
{
	if (nstgs_ == 0)
	{
		std::cerr << "No stochastic problem is loaded." << std::endl;
		return DSP_RTN_ERR;
	}
	if (ncols_[0] == 0)
	{
		std::cerr << "No stochastic programming variable is loaded." << std::endl;
		return DSP_RTN_ERR;
	}
	if (nscen_ == 0)
	{
		std::cerr << "No scenario is loaded." << std::endl;
		return DSP_RTN_ERR;
	}

	wass_eps_ = pow(eps, lp_norm);
	isdro_ = true;

	/** Count the number of reference scenarios.
	 * We take all scenarios with nonzero probabilities as the references (a.k.a. 
	 * empirical observations).
	 */
	nrefs_ = 0;
	for (int s = 0; s < nscen_; ++s)
	{
		if (prob_[s] > 0)
			nrefs_++;
	}

	refs_probability_ = new double[nrefs_];
	wass_dist_ = new double *[nrefs_];

	/** Compute the Wasserstein distances of lp_norm,
	 * and assign the reference probabilities.
	 */
	for (int s = 0, r = 0; s < nscen_; ++s)
	{
		if (prob_[s] > 0)
		{
			refs_probability_[r] = prob_[s];
			wass_dist_[r] = new double[nscen_];
			CoinZeroN(wass_dist_[r], nscen_);
			assert(mat_scen_[s]->getNumRows() == nrows_[1]);

			/** Compute the distances.
			 * 
			 * TODO: Can we do in parallel?
			 * The relevant issues need addressed first: 
			 * - https://github.com/kibaekkim/DSPopt.jl/issues/14
			 * - https://github.com/Argonne-National-Laboratory/DSP/issues/50
			 */
			for (int ss = 0; ss < nscen_; ++ss)
			{
				for (int j = 0; j < ncols_[1]; ++j)
				{
					wass_dist_[r][ss] += pow(fabs((*obj_scen_[s])[j] - (*obj_scen_[ss])[j]), 2);
					if ((*clbd_scen_[s])[j] > -1.e+20 && (*clbd_scen_[ss])[j] > -1.e+20)
						wass_dist_[r][ss] += pow(fabs((*clbd_scen_[s])[j] - (*clbd_scen_[ss])[j]), 2);
					if ((*cubd_scen_[s])[j] < 1.e+20 && (*cubd_scen_[ss])[j] < 1.e+20)
						wass_dist_[r][ss] += pow(fabs((*cubd_scen_[s])[j] - (*cubd_scen_[ss])[j]), 2);
				}
				for (int i = 0; i < nrows_[1]; ++i)
				{
					if ((*rlbd_scen_[s])[i] > -1.e+20 && (*rlbd_scen_[ss])[i] > -1.e+20)
						wass_dist_[r][ss] += pow(fabs((*rlbd_scen_[s])[i] - (*rlbd_scen_[ss])[i]), 2);
					if ((*rubd_scen_[s])[i] < 1.e+20 && (*rubd_scen_[ss])[i] < 1.e+20)
						wass_dist_[r][ss] += pow(fabs((*rubd_scen_[s])[i] - (*rubd_scen_[ss])[i]), 2);
					for (int j = 0; j < mat_scen_[s]->getNumCols(); ++j)
					{
						wass_dist_[r][ss] += pow(fabs(mat_scen_[s]->getCoefficient(i, j) - mat_scen_[ss]->getCoefficient(i, j)), 2);
					}
				}
				wass_dist_[r][ss] = pow(wass_dist_[r][ss], lp_norm / 2.0);
			}
			r++;
		}
	}

	/** scaling vector */
	double scaling_constant = pow(wass_eps_, 2);
	for (int r = 0; r < nrefs_; ++r)
		for (int s = 0; s < nscen_; ++s)
			scaling_constant += pow(wass_dist_[r][s], 2);
	scaling_constant = sqrt(scaling_constant);
	wass_eps_ /= scaling_constant;
	for (int r = 0; r < nrefs_; ++r)
		for (int s = 0; s < nscen_; ++s)
			wass_dist_[r][s] /= scaling_constant;

	/** Quadratic equations
	 * TODO: The quadratic objective function and constraints need to be considered.
	 * - https://github.com/Argonne-National-Laboratory/DSP/issues/155
	 */

	printf("[DRO] Set %d reference scenarios.\n", nrefs_);
	printf("[DRO] Computed the Wasserstein distances with %f-norm.\n", lp_norm);

	return DSP_RTN_OK;
}

/** split core matrix row for a given stage */
CoinPackedVector * StoModel::splitCoreRowVec(
		int i,  /**< row index */
		int stg /**< stage */)
{
	int j, ind;

	/** create empty vector */
	CoinPackedVector * split = new CoinPackedVector;

	/** reserve memory */
	split->reserve(ncols_[stg]);

	/** insert elements */
	for (j = 0; j < rows_core_[i]->getNumElements(); ++j)
	{
		ind = rows_core_[i]->getIndices()[j] - cstart_[stg];
		if (ind >= 0 && ind < ncols_[stg])
		{
			split->insert(ind, rows_core_[i]->getElements()[j]);
		}
	}

	return split;
}

/** copy core column lower bounds */
void StoModel::copyCoreColLower(double * clbd, int stg)
{
	CoinCopyN(clbd_core_[stg], ncols_[stg], clbd);
}

/** copy core column upper bounds */
void StoModel::copyCoreColUpper(double * cubd, int stg)
{
	CoinCopyN(cubd_core_[stg], ncols_[stg], cubd);
}

/** copy core objective coefficients */
void StoModel::copyCoreObjective(double * obj, int stg)
{
	for (int j = 0; j < ncols_[stg]; ++j)
	{
		obj[j] = obj_core_[stg][j];
	}
}

void StoModel::copyCoreQuadrativeObjective(CoinPackedMatrix * &qobj_coupling, CoinPackedMatrix * &qobj_ncoupling, int stg)
{
	vector<int> colidx;
	vector<int> rowidx;
	vector<double> elements;
	int numq=qobj_core_[stg]->getNumElements();
	const CoinBigIndex * start = qobj_core_[stg]->getVectorStarts();
	vector<int> colidx_coupling;
	vector<int> rowidx_coupling;
	vector<double> elements_coupling;
	int numq_coupling = 0;
	vector<int> colidx_ncoupling;
	vector<int> rowidx_ncoupling;
	vector<double> elements_ncoupling;
	int numq_ncoupling = 0;

	int rowcount=0;
	for (int i=0; i<qobj_core_[stg]->getMajorDim();i++){
		if (start[i+1]-start[i]>0){
				for (int k=0; k<start[i+1]-start[i]; k++){
					rowidx.push_back(rowcount);
				}
		}
		rowcount++;
	}

	for (int i=0; i<numq; i++){
		colidx.push_back(qobj_core_[stg]->getIndices()[i]);
		elements.push_back(qobj_core_[stg]->getElements()[i]);
	}

	// security check
	assert(numq==rowidx.size());
	assert(numq==colidx.size());
	assert(numq==elements.size());
	for (int i=0; i<numq; i++){
		if (rowidx[i]>=ncols_[stg-1] && colidx[i]>=ncols_[stg-1]){
			rowidx_ncoupling.push_back(rowidx[i]-ncols_[stg-1]);
			colidx_ncoupling.push_back(colidx[i]-ncols_[stg-1]);
			elements_ncoupling.push_back(elements[i]);
			numq_ncoupling++;
		}
		else{
			rowidx_coupling.push_back(rowidx[i]);
			colidx_coupling.push_back(colidx[i]);
			elements_coupling.push_back(elements[i]);
			numq_coupling++;
		}
	}

	if(numq_ncoupling!=0){
		qobj_ncoupling = new CoinPackedMatrix(false, &rowidx_ncoupling[0], &colidx_ncoupling[0], &elements_ncoupling[0], numq_ncoupling);
	}
	else{
		qobj_ncoupling = NULL;
	}
	
	if(numq_coupling!=0){
		qobj_coupling = new CoinPackedMatrix(false, &rowidx_coupling[0], &colidx_coupling[0], &elements_coupling[0], numq_coupling);
	}
	else{
		qobj_coupling = NULL;
	}
}

/** copy core column types */
void StoModel::copyCoreColType(char * ctype, int stg)
{
	CoinCopyN(ctype_core_[stg], ncols_[stg], ctype);
}

/** copy core row lower bounds */
void StoModel::copyCoreRowLower(double * rlbd, int stg)
{
	CoinCopyN(rlbd_core_[stg], nrows_[stg], rlbd);
}

/** copy core row upper bounds */
void StoModel::copyCoreRowUpper(double * rubd, int stg)
{
	CoinCopyN(rubd_core_[stg], nrows_[stg], rubd);
}

/** combine random matrix row for a given scenario */
void StoModel::combineRandRowVec(
		CoinPackedVector * row, /**< core row vector */
		int i,                  /**< row index */
		int scen                /**< scenario index */)
{
	int j;
	int pos = 0; /** position to sparse row vector elements */

	/** random matrix row info */
	CoinBigIndex   start = mat_scen_[scen]->getVectorStarts()[i];
	const int *    ind   = mat_scen_[scen]->getIndices() + start;
	const double * elem  = mat_scen_[scen]->getElements() + start;

	for (j = 0; j < mat_scen_[scen]->getVectorSize(i); ++j)
	{
		while (row->getIndices()[pos] < ind[j] && pos < row->getNumElements())
		{
			pos++;
		}
		if (row->getIndices()[pos] == ind[j])
		{
			row->setElement(pos, elem[j]);
		}
	}
}

/** combine random matrix row for a given stage and scenario */
void StoModel::combineRandRowVec(
		CoinPackedVector * row, /**< core row vector */
		int i,                  /**< row index */
		int stg,                /**< stage index */
		int scen                /**< scenario index */)
{
	int j, col;
	int pos = 0; /** position to sparse row vector elements */

	/** random matrix row info */
	CoinBigIndex   start = mat_scen_[scen]->getVectorStarts()[i];
	const int *    ind   = mat_scen_[scen]->getIndices() + start;
	const double * elem  = mat_scen_[scen]->getElements() + start;

	for (j = 0; j < mat_scen_[scen]->getVectorSize(i); ++j)
	{
		col = ind[j] - cstart_[stg];
		if (col >= 0 && col < ncols_[stg])
		{
			while (row->getIndices()[pos] < col && pos < row->getNumElements())
			{
				pos++;
			}
			if (row->getIndices()[pos] == col)
			{
				row->setElement(pos, elem[j]);
			}
		}
	}
}	

/** combine random column lower bounds */
void StoModel::combineRandColLower(double * clbd, int stg, int scen)
{
	for (int i = 0; i < clbd_scen_[scen]->getNumElements(); ++i)
	{
		int j = clbd_scen_[scen]->getIndices()[i] - cstart_[stg];
		if (j >= 0 && j < ncols_[stg])
		{
			clbd[j] = clbd_scen_[scen]->getElements()[i];
		}
	}
}

/** combine random column upper bounds */
void StoModel::combineRandColUpper(double * cubd, int stg, int scen)
{
	for (int i = 0; i < cubd_scen_[scen]->getNumElements(); ++i)
	{
		int j = cubd_scen_[scen]->getIndices()[i] - cstart_[stg];
		if (j >= 0 && j < ncols_[stg])
		{
			cubd[j] = cubd_scen_[scen]->getElements()[i];
		}
	}
}

/** combine random objective coefficients */
void StoModel::combineRandObjective(double * obj, int stg, int scen, bool adjustProbability)
{
	int i, j;
	for (i = 0; i < obj_scen_[scen]->getNumElements(); ++i)
	{
		j = obj_scen_[scen]->getIndices()[i] - cstart_[stg];
		if (j >= 0 && j < ncols_[stg])
		{
			obj[j] = obj_scen_[scen]->getElements()[i];
		}
	}

	if (adjustProbability)
	{ 
		for (j = ncols_[stg] - 1; j >= 0; --j)
		{
			obj[j] *= prob_[scen];
		}
	}
}

void StoModel::combineRandQuadraticObjective(CoinPackedMatrix * &qobj_coupling, CoinPackedMatrix * &qobj_ncoupling, int stg, int scen, bool adjustProbability)
{
	int i, j, s;
	int numelements = qobj_scen_[scen]->getNumElements();
	const int * colindex=qobj_scen_[scen]->getIndices();
	const CoinBigIndex * start = qobj_scen_[scen]->getVectorStarts();
	const int dim=qobj_scen_[scen]->getMajorDim();
	const double * elements=qobj_scen_[scen]->getElements();

	int * rowindex;
	rowindex = new int [numelements];
	
	std::vector<int> newrowindex_coupling, newcolindex_coupling;
	std::vector<double> newelements_coupling;

	std::vector<int> newrowindex_ncoupling, newcolindex_ncoupling;
	std::vector<double> newelements_ncoupling;

	
	int idx=0;
	/** convert start into row index */
	for (int k=0; k<dim; k++){
		for (i=0; i<start[k+1]-start[k];i++){
			rowindex[idx]=k;
			idx++;
		}
	}

	int sidx=stg-1;
	if (sidx<0) sidx=0;
	for (i = 0; i < numelements; ++i)
	{
		if (colindex[i]>=ncols_[sidx] && rowindex[i]>=ncols_[sidx]){
			newrowindex_ncoupling.push_back(rowindex[i]-ncols_[sidx]);
			newcolindex_ncoupling.push_back(colindex[i]-ncols_[sidx]);
			if (adjustProbability){
				newelements_ncoupling.push_back(elements[i]*prob_[scen]);
			}
			else{
				newelements_ncoupling.push_back(elements[i]);
			}
		}
		else{
			newrowindex_coupling.push_back(rowindex[i]);
			newcolindex_coupling.push_back(colindex[i]);
			if (adjustProbability){
				newelements_coupling.push_back(elements[i]*prob_[scen]);
			}
			else{
				newelements_coupling.push_back(elements[i]);
			}
		}
	}
	//printf("newcolindex_coupling.size() = %d\n", newcolindex_coupling.size());
	//for (int m=0; m<newcolindex_ncoupling.size(); m++){
	//	printf("newcolindex_ncoupling[%d]=%d\n", m, newcolindex_ncoupling[m]);
	//}
	if (newelements_coupling.size()==0){
		qobj_coupling=NULL;
	}
	else{
		qobj_coupling=new CoinPackedMatrix(false, &newcolindex_coupling[0], &newrowindex_coupling[0], &newelements_coupling[0], newelements_coupling.size());
	}
	if (newelements_ncoupling.size()==0){
		qobj_coupling=NULL;
	}
	else{
		qobj_ncoupling=new CoinPackedMatrix(false, &newcolindex_ncoupling[0], &newrowindex_ncoupling[0], &newelements_ncoupling[0], newelements_ncoupling.size());
	}
	//qobj_ncoupling=new CoinPackedMatrix(false, &newcolindex_ncoupling[0], &newrowindex_ncoupling[0], &newelements_ncoupling[0], newelements_ncoupling.size());
	//PRINT_ARRAY_MSG(newcolindex_ncoupling.size(), &newcolindex_ncoupling[0], "elements in newcolindex_ncoupling");
	//PRINT_ARRAY_MSG(qobj_ncoupling->getNumElements(), qobj_ncoupling->getElements(), "elements in qobj_ncoupling");

	delete[] rowindex;
}

/** combine random row lower bounds */
void StoModel::combineRandRowLower(double * rlbd, int stg, int scen)
{
	for (int i = 0; i < rlbd_scen_[scen]->getNumElements(); ++i)
	{
		int j = rlbd_scen_[scen]->getIndices()[i] - rstart_[stg];
		if (j >= 0 && j < nrows_[stg])
		{
			rlbd[j] = rlbd_scen_[scen]->getElements()[i];
		}
	}
}

/** combine random row upper bounds */
void StoModel::combineRandRowUpper(double * rubd, int stg, int scen)
{
	for (int i = 0; i < rubd_scen_[scen]->getNumElements(); ++i)
	{
		int j = rubd_scen_[scen]->getIndices()[i] - rstart_[stg];
		if (j >= 0 && j < nrows_[stg])
		{
			rubd[j] = rubd_scen_[scen]->getElements()[i];
		}
	}
}

/** shift vector indices by offset */
void StoModel::shiftVecIndices(
		CoinPackedVector * vec, /**< vector to shift indices */
		int offset,             /**< offset by which indices are shifted */
		int start               /**< index only after which indices are shifted */)
{
	shiftVecIndices(vec->getNumElements(), vec->getIndices(), offset, start);
}

/** shift vector indices by offset */
void StoModel::shiftVecIndices(
		int size,           /**< size of vecind */
		int * vecind, /**< vector indices to shift */
		int offset,         /**< offset by which indices are shifted */
		int start           /**< index only after which indices are shifted */)
{
	if (offset == 0) return;
	for (int i = size - 1; i >= 0; --i)
	{
		if (vecind[i] >= start)
		{
			vecind[i] += offset;
		}
	}
}

#if 0
/** add branching object */
void StoModel::addBranchingHyperplane(int nzcnt, int * indices, double * values, int priority)
{
	CoinPackedVector * vec = new CoinPackedVector;
	vec->setVector(nzcnt, indices, values);
	branchingHyperplanes_.push_back(
			make_pair(vec, priority));
	//printf("number of branching hyperplanes %d\n", (int)branchingHyperplanes_.size());
}
#endif

void StoModel::__printData()
{
	char tmpstr[128];

	printf("\n### BEGINNING of printing StoModel data ###\n\n");

	printf("nscen_ %d\n", nscen_);
	printf("nstgs_ %d\n", nstgs_);
	PRINT_ARRAY_MSG(nstgs_, nrows_, "nrows_")
	PRINT_ARRAY_MSG(nstgs_, ncols_, "ncols_")
	PRINT_ARRAY_MSG(nstgs_, rstart_, "rstart_")
	PRINT_ARRAY_MSG(nstgs_, cstart_, "cstart_")

	printf("\n### Core data ###\n\n");
	printf("nrows_core_ %d\n", nrows_core_);
	printf("ncols_core_ %d\n", ncols_core_);
#if 1
	for (int i = 0; i < nrows_core_; ++i)
	{
		sprintf(tmpstr, "Core matrix row %d", i);
		PRINT_COIN_PACKED_VECTOR_MSG((*rows_core_[i]), tmpstr)
	}
#endif
	for (int t = 0; t < nstgs_; ++t)
	{
		printf("\n### Stage %d core data ###\n\n", t);
		PRINT_ARRAY_MSG(ncols_[t], clbd_core_[t], "clbd_core_")
		PRINT_ARRAY_MSG(ncols_[t], cubd_core_[t], "cubd_core_")
		PRINT_ARRAY_MSG(ncols_[t], obj_core_[t], "obj_core_")
		PRINT_ARRAY_MSG(ncols_[t], ctype_core_[t], "ctype_core_")
		PRINT_ARRAY_MSG(nrows_[t], rlbd_core_[t], "rlbd_core_")
		PRINT_ARRAY_MSG(nrows_[t], rubd_core_[t], "rubd_core_")
		printf("\n### quadratic objective coefficient ###\n\n");
		if (qobj_core_[t])
		{
			printf("isColOrdered %d\n", qobj_core_[t]->isColOrdered());
			PRINT_ARRAY_MSG(qobj_core_[t]->getMajorDim(), qobj_core_[t]->getVectorStarts(), "VectorStarts")
			PRINT_SPARSE_ARRAY_MSG(
					qobj_core_[t]->getNumElements(),
					qobj_core_[t]->getIndices(),
					qobj_core_[t]->getElements(),
					"Elements")
		}
	}

	for (int s = 0; s < nscen_; ++s)
	{
		printf("\n### Scenario %d data ###\n\n", s);
		printf("probability %E\n", prob_[s]);
		printf("stage map %d\n", scen2stg_[s]);
#if 1
		printf("=== BEGINNING of CoinPackedMatrix mat_scen_[%d] ===\n", s);
		if (mat_scen_[s])
		{
			printf("isColOrdered %d\n", mat_scen_[s]->isColOrdered());
			PRINT_ARRAY_MSG(mat_scen_[s]->getMajorDim(), mat_scen_[s]->getVectorStarts(), "VectorStarts")
			PRINT_SPARSE_ARRAY_MSG(
					mat_scen_[s]->getNumElements(),
					mat_scen_[s]->getIndices(),
					mat_scen_[s]->getElements(),
					"Elements")
		}
		printf("=== END of CoinPackedMatrix mat_scen_[%d] ===\n", s);

		printf("=== BEGINNING of CoinPackedMatrix qobj_scen_[%d] ===\n", s);
		if (qobj_scen_[s])
		{
			printf("isColOrdered %d\n", qobj_scen_[s]->isColOrdered());
			PRINT_ARRAY_MSG(qobj_scen_[s]->getMajorDim(), qobj_scen_[s]->getVectorStarts(), "VectorStarts")
			PRINT_SPARSE_ARRAY_MSG(
					qobj_scen_[s]->getNumElements(),
					qobj_scen_[s]->getIndices(),
					qobj_scen_[s]->getElements(),
					"Elements")
		}
		printf("=== END of CoinPackedMatrix qobj_scen_[%d] ===\n", s);
#endif
//		PRINT_COIN_PACKED_VECTOR_MSG((*clbd_scen_[s]), "clbd_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*cubd_scen_[s]), "cubd_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*obj_scen_[s]), "obj_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*rlbd_scen_[s]), "rlbd_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*rubd_scen_[s]), "rubd_scen_")
	}

	printf("\n### END of printing StoModel data ###\n");
}

DSP_RTN_CODE StoModel::printQuadRows(const int s)
{
	BGN_TRY_CATCH

	if (s < 0) {
		for (int i = 0; i < qc_row_core_->nqrows; i++) 
		{
			cout << "Core " << i << "th quad constr: ";
			for (int lt = 0; lt < qc_row_core_->linnzcnt[i]; lt++)
			{
				cout << qc_row_core_->linval[i][lt] << " x" << qc_row_core_->linind[i][lt] << " + ";
			}
			for (int qt = 0; qt < qc_row_core_->quadnzcnt[i]-1; qt++)
			{
				cout << qc_row_core_->quadval[i][qt] << " x" << qc_row_core_->quadrow[i][qt] << " x" << qc_row_core_->quadcol[i][qt] << " + ";
			}
			cout << qc_row_core_->quadval[i][qc_row_core_->quadnzcnt[i]-1] << " x" << qc_row_core_->quadrow[i][qc_row_core_->quadnzcnt[i]-1] << " x" << qc_row_core_->quadcol[i][qc_row_core_->quadnzcnt[i]-1];
			if (qc_row_core_->sense[i] == 'L')
				cout << " <= " << qc_row_core_->rhs[i] << endl;
			else 
				cout << " >= " << qc_row_core_->rhs[i] << endl;
		}
	} else {
		for (int i = 0; i < qc_row_scen_[s].nqrows; i++) 
		{
			cout << "Scen " << s << "th " << i << "th quad constr: ";
			for (int lt = 0; lt < qc_row_scen_[s].linnzcnt[i]; lt++)
			{
				cout << qc_row_scen_[s].linval[i][lt] << " x" << qc_row_scen_[s].linind[i][lt] << " + ";
			}
			for (int qt = 0; qt < qc_row_scen_[s].quadnzcnt[i]-1; qt++)
			{
				cout << qc_row_scen_[s].quadval[i][qt] << " x" << qc_row_scen_[s].quadrow[i][qt] << " x" << qc_row_scen_[s].quadcol[i][qt] << " + ";
			}
			cout << qc_row_scen_[s].quadval[i][qc_row_scen_[s].quadnzcnt[i]-1] << " x" << qc_row_scen_[s].quadrow[i][qc_row_scen_[s].quadnzcnt[i]-1] << " x" << qc_row_scen_[s].quadcol[i][qc_row_scen_[s].quadnzcnt[i]-1];
			if (qc_row_scen_[s].sense[i] == 'L')
				cout << " <= " << qc_row_scen_[s].rhs[i] << endl;
			else 
				cout << " >= " << qc_row_scen_[s].rhs[i] << endl;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE StoModel::printQuadRows(const QuadRowData *qc)
{
	BGN_TRY_CATCH

	for (int i = 0; i < qc->nqrows; i++) 
	{
		cout << i << "th quad constr: ";

		for (int lt = 0; lt < qc->linnzcnt[i]; lt++)
		{
			cout << qc->linval[i][lt] << " x" << qc->linind[i][lt] << " + ";
		}
		for (int qt = 0; qt < qc->quadnzcnt[i]-1; qt++)
		{
			cout << qc->quadval[i][qt] << " x" << qc->quadrow[i][qt] << " x" << qc->quadcol[i][qt] << " + ";
		}
		cout << qc->quadval[i][qc->quadnzcnt[i]-1] << " x" << qc->quadrow[i][qc->quadnzcnt[i]-1] << " x" << qc->quadcol[i][qc->quadnzcnt[i]-1];
		if (qc->sense[i] == 'L')
			cout << " <= " << qc->rhs[i] << endl;
		else 
			cout << " >= " << qc->rhs[i] << endl;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

double StoModel::getWassersteinDist(int i, int j) {
	if (!isdro_) return 0.0;
	if (i >= nrefs_) {
		printf("Reference index (%d) is out of range (%d).\n", i, nrefs_);
		return 0.0;
	}
	if (j >= nscen_) {
		printf("Scenario index (%d) is out of range (%d).\n", j, nscen_);
		return 0.0;
	}
	return wass_dist_[i][j];
}

double StoModel::getReferenceProbability(int i) {
	if (refs_probability_ == NULL) return 0.0;
	if (!isdro_) return 0.0;
	if (i >= nrefs_) {
		printf("Reference index (%d) is out of range (%d).\n", i, nrefs_);
		return 0.0;
	}
	return refs_probability_[i];
}
