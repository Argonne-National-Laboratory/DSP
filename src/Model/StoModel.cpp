/*
 * StoModel.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: kibaekkim
 */
// #define DSP_DEBUG
#include "CoinHelperFunctions.hpp"
#include "Utility/DspMessage.h"
#include "StoModel.h"
#include <sstream>

extern "C" void dpotrf_(char *uplo,int *jb, double *A, int *lda, int *info); 

/* Prints matrix in column-major format. */
static void show_matrix(const double* A, const int n)
{
    int i = 0, j = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%2.5f ", A[i * n + j]);
        }
        printf("\n");
    }
}

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
	
	if (rhs.qc_row_scen_ != NULL) 
	{
		qc_row_scen_ = new QuadRowData * [nscen_];
	} 
	else 
	{
		qc_row_scen_ = NULL;
	}

	for (int i = nscen_ - 1; i >= 0; --i)
	{
		if (qc_row_scen_) {
			qc_row_scen_[i] = rhs.qc_row_scen_[i];
		}
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
	FREE_2D_PTR(nscen_, qc_row_scen_);
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
			int idx = core->getBinaryIndices()[j] - cstart_[i];
			ctype_core_[i][idx] = 'B';
			if (clbd_core_[i][idx] < 0.0)
				clbd_core_[i][idx] = 0.0;
			if (clbd_core_[i][idx] > 1.0)
				cubd_core_[i][idx] = 1.0;
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
		mat_scen_[i]  = new CoinPackedMatrix(false, ncols_core_, nrows_[stg], node->getNumMatrixElements(), 
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

void StoModel::setProbability(double *probability)
{
	CoinCopyN(probability, nscen_, prob_);
}

/* split a string */
template <class Container>
void split(const std::string& str, Container& cont)
{
	cont.clear();
    istringstream iss(str);
    copy(istream_iterator<std::string>(iss),
         istream_iterator<std::string>(),
         back_inserter(cont));
}

void addAffineRowOfL(double * L, int n, vector<int> indices, CoinPackedVector ** &rows_core_temp, int cstart_of_col_to_add, int rstart_of_row_to_add, int linnzcnt, int *linind, double *linval)
{
	int nzcnt;
	int j, k;

	/* constraint for jth row of Q[i] */
	for (j = 0; j < n; j++) 
	{
		nzcnt = 0;
		for (k = 0; k < n; k++) 
		{
			if (fabs(L[n * j + k]) > 1e-8) {
				nzcnt++;
			}
		}

		/* term for additional variables */
		nzcnt += 1;
		int * idx = new int [nzcnt];
		double * val = new double [nzcnt];
		int pos = 0;
		for (k = 0; k < n; k++) 
		{
			if (fabs(L[n * j + k]) > 1e-8) {
				idx[pos] = indices[k];
				val[pos] = L[n * j + k];
				pos++;
			}
		}

		assert(pos == nzcnt - 1);
		idx[nzcnt-1] = cstart_of_col_to_add + j;
		val[nzcnt-1] = -1;
		rows_core_temp[rstart_of_row_to_add + j] = new CoinPackedVector(nzcnt, idx, val);	

		FREE_ARRAY_PTR(idx);
		FREE_ARRAY_PTR(val);
	}

	/* last two additional rows */
	int * ind = new int [linnzcnt+1];
	double * val = new double [linnzcnt+1];
	double * val_n = new double [linnzcnt+1];
	for (j = 0; j < linnzcnt; j++) 
	{
		ind[j] = linind[j];
		val[j] = linval[j];
		val_n[j] = -linval[j];
	}
	ind[linnzcnt] = cstart_of_col_to_add + n;
	val[linnzcnt] = -2;
	val_n[linnzcnt] = -2;

	int row_idx = rstart_of_row_to_add + n;
	
	rows_core_temp[row_idx] = new CoinPackedVector(linnzcnt+1, ind, val);

	row_idx = rstart_of_row_to_add + n + 1;
	ind[linnzcnt] = cstart_of_col_to_add + n + 1;
	rows_core_temp[row_idx] = new CoinPackedVector(linnzcnt+1, ind, val_n);

	FREE_ARRAY_PTR(ind);
	FREE_ARRAY_PTR(val);
	FREE_ARRAY_PTR(val_n);
}
void StoModel::getL(double * &Q, int quadnzcnt, int *quadcol, int *quadrow, double *quadval, vector<int> &indices, int &n)
{
	int idx;
	int i, j, k;

	vector<int>::iterator itr;

	for (j = 0; j < quadnzcnt; j++) 
	{	
		for (k = 0; k < 2; k++)
		{
			if (k == 0)
				idx = quadcol[j];
			else 
				idx = quadrow[j];

			itr = find (indices.begin(), indices.end(), idx);
			if (itr == indices.end())
				indices.push_back(idx);
		}
	}
	sort(indices.begin(), indices.end()); 

	n = indices.size();

	/* make symmetric QCMATIX */
	int col, row;
	Q = new double [n * n] {0};
	for (j = 0; j < quadnzcnt; j++) 
	{
		itr = find (indices.begin(), indices.end(), quadcol[j]);
		col = itr - indices.begin();
		itr = find (indices.begin(), indices.end(), quadrow[j]);
		row = itr - indices.begin();
		/* make it symmetric */
		Q[n * row + col] += 0.5 * quadval[j];
		Q[n * col + row] += 0.5 * quadval[j];
	}
#ifdef DSP_DEBUG
show_matrix(Q, n);
#endif
	/* remove zero columns or rows */
	vector<int> indices_temp;
	int n_temp;

	for (i = 0; i < n; i++)
	{
		bool is_null_row = true;
		for (j = 0; j < n; j++) 
		{
			if (fabs(Q[n * i + j]) > 1e-6)
			{
				is_null_row = false;
				break;
			}
		}
		if (!is_null_row) indices_temp.push_back(i);
	}
	if (indices_temp.size() < n) 
	{
		n_temp = indices_temp.size();
		double * Q_temp = NULL;
		
		if (n_temp > 0)
		{
			Q_temp = new double [n_temp * n_temp] {0};

			for (i = 0; i < n_temp; i++)
			{
				for (j = 0; j < n_temp; j++)
				{
					Q_temp[i * n_temp + j] = Q[indices_temp[i] * n + indices_temp[j]];
					Q_temp[j * n_temp + i] = Q[indices_temp[j] * n + indices_temp[i]];
				}
			}
#ifdef DSP_DEBUG
show_matrix(Q_temp, n_temp);
#endif
			/* Cholesky Decomposition */
			int lda = n_temp;     /* leading dimension of array */
			char uplo = 'L'; /* the lower diagonal matrix L */
			int info = 0;    /* returns non-zero value for error */
			dpotrf_(&uplo, &n_temp, Q_temp, &lda, &info);
			// int * piv = new int [n_temp];
			// int rank;
			// double tol = 1e-8;
			// double * work = new double [2*n_temp];
			// dpstrf_(&uplo, &n_temp, Q_temp, &lda, &info, piv, &tol, work, &info);
			// assert(info == 0);
		}
		
		delete Q;
		Q = new double [n * n] {0};
		
		for (i = 0; i < n_temp; i++)
		{
			for (j = 0; j < n_temp; j++)
			{
				Q[indices_temp[i] * n + indices_temp[j]] = Q_temp[i * n_temp + j];
				Q[indices_temp[j] * n + indices_temp[i]] = Q_temp[j * n_temp + i];
			}
		}
	} 
	else 
	{
		/* Cholesky Decomposition */
		int lda = n;     /* leading dimension of array */
		char uplo = 'L'; /* the lower diagonal matrix L */
		int info = 0;    /* returns non-zero value for error */
		dpotrf_(&uplo, &n, Q, &lda, &info);
		// int * piv = new int [n];
		// int rank;
		// double tol = 1e-8;
		// double * work = new double [2*n];
		// dpstrf_(&uplo, &n, Q, &lda, &info, piv, &tol, work, &info);
		// assert(info == 0);
	}
	/* Note that because m1 is over-written by LAPACK, you need to 0 out the 
	lower diagonal entries yourself, since LAPACK will not do it for you. */
	for (i = 0; i < n; i++)
    	for (j = 0; j < n; j++)
        	if (i < j)
				Q[j * n + i] = 0.0;
#ifdef DSP_DEBUG
show_matrix(Q, n);
#endif
}
DSP_RTN_CODE StoModel::chgToSocp(vector<int> &qc_rstart)
{
	BGN_TRY_CATCH

	int i, j, k, s;

	/* nqrows: total qc core rows 
	 * qcrows_stg: number of qc rows in each stages */
	int nqrows = qc_row_core_->nqrows;
	int * nqrows_stg = new int [nstgs_];
	for (i = 0; i < nstgs_; i++)
	{
		if (i < nstgs_ - 1)
			nqrows_stg[i] = qc_rstart[i+1] - qc_rstart[i];
		else 
			nqrows_stg[i] = nqrows - qc_rstart[i];
	}
	/* allocate memory for the lower-diagonal matrices of QcMatrix 
	 * L is a n by n matrix for some n, indices: n number of indices for L */ 
	double ** L_core = new double * [nqrows];
	int * n_core = new int [nqrows] {0};
	vector<vector<int>> indices_core(nqrows);
	
	double *** L_scen = new double ** [nscen_];
	int ** n_scen = new int * [nscen_];
	vector<vector<vector<int>>> indices_scen(nscen_, vector<vector<int>>(nqrows_stg[1]));
cout << "reach 735" << endl;
	for (s = 0; s < nscen_; s++) 
	{
		L_scen[s] = new double * [nqrows_stg[1]];
		n_scen[s] = new int [nqrows_stg[1]] {0};
	}

	/* for each QcMatrix, construct the lower diagonal matrix L */
	int naddrows = 0;
	int naddcols = 0;
	for (i = 0; i < nqrows; i++) 
	{
		/* for core qc rows */
		int linnzcnt = qc_row_core_->linnzcnt[i];
		int quadnzcnt = qc_row_core_->quadnzcnt[i];
		
		int * quadcol = qc_row_core_->quadcol[i];
		int * quadrow = qc_row_core_->quadrow[i];
		double * quadval = qc_row_core_->quadval[i];

		getL(L_core[i], quadnzcnt, quadcol, quadrow, quadval, indices_core[i], n_core[i]);

		/* for scenario qc rows */
		if (i >= qc_rstart[1]) 
		{
			int row_idx = i - qc_rstart[1];
			
			for (int s = 0; s < nscen_; s++)
			{
#ifdef DSP_DEBUG
cout << "L_scen[" << s << "]: " << endl; 
#endif
				linnzcnt = qc_row_scen_[s]->linnzcnt[row_idx];
				quadnzcnt = qc_row_scen_[s]->quadnzcnt[row_idx];
				
				quadcol = qc_row_scen_[s]->quadcol[row_idx];
				quadrow = qc_row_scen_[s]->quadrow[row_idx];
				quadval = qc_row_scen_[s]->quadval[row_idx];

				getL(L_scen[s][row_idx], quadnzcnt, quadcol, quadrow, quadval, indices_scen[s][row_idx], n_scen[s][row_idx]);

				if (n_core[i] != n_scen[s][row_idx])
				{
					cerr << "All terms for quadratic constraints should be declared in the core file." << endl;
					return DSP_RTN_ERR;
				}
			}
		}

		/* number of added linear constraints and variables to the core */
		naddrows += n_core[i] + 2;
		naddcols += n_core[i] + 2;
	}

	/* resize and modify smps data */
	int nrows_core_temp = nrows_core_ + naddrows;
	int ncols_core_temp = ncols_core_ + naddcols;

	int * nrows_temp = new int [nstgs_];
	int * ncols_temp  = new int [nstgs_];
	int * rstart_temp  = new int [nstgs_];
	int * cstart_temp  = new int [nstgs_];

	/* copy original smps data */
	CoinCopyN(nrows_, nstgs_, nrows_temp);
	CoinCopyN(ncols_, nstgs_, ncols_temp);
	CoinCopyN(rstart_, nstgs_, rstart_temp);
	CoinCopyN(cstart_, nstgs_, cstart_temp);
	
	/* reflect the added linear rows and cols */
	for (i = 0; i < qc_rstart[1]; i++)
	{
		nrows_temp[0] += n_core[i] + 2;
		ncols_temp[0] += n_core[i] + 2;

		rstart_temp[1] += n_core[i] + 2;
		cstart_temp[1] += n_core[i] + 2;
	}
	for (i = qc_rstart[1]; i < nqrows; i++)
	{
		nrows_temp[1] += n_core[i] + 2;
		ncols_temp[1] += n_core[i] + 2;
	}

	assert(nrows_core_temp == nrows_temp[0] + nrows_temp[1]);
	assert(ncols_core_temp == ncols_temp[0] + ncols_temp[1]);

#ifdef DSP_DEBUG
for (j = 0; j < nstgs_; j++)
{
	cout << "nrows_[" << j << "]: " << nrows_[j] << endl;
	cout << "nrows_temp[" << j << "]: " << nrows_temp[j] << endl;
}
for (j = 0; j < nstgs_; j++)
{
	cout << "ncols_[" << j << "]: " << ncols_[j] << endl;
	cout << "ncols_temp[" << j << "]: " << ncols_temp[j] << endl;
}
for (j = 0; j < nstgs_; j++)
{
	cout << "rstart_[" << j << "]: " << rstart_[j] << endl;
	cout << "rstart_temp[" << j << "]: " << rstart_temp[j] << endl;
}
for (j = 0; j < nstgs_; j++)
{
	cout << "cstart_[" << j << "]: " << cstart_[j] << endl;
	cout << "cstart_temp[" << j << "]: " << cstart_temp[j] << endl;
}
#endif
	
	double ** obj_core_temp = new double * [nstgs_];
	double ** clbd_core_temp = new double * [nstgs_];
	double ** cubd_core_temp = new double * [nstgs_];
	char ** ctype_core_temp = new char * [nstgs_];
	double ** rlbd_core_temp = new double * [nstgs_];
	double ** rubd_core_temp = new double * [nstgs_];

	for (j = 0; j < nstgs_; j++)
	{
		obj_core_temp[j] = new double [ncols_temp[j]] {0};
		CoinCopyN(obj_core_[j], ncols_[j], obj_core_temp[j]);

		clbd_core_temp[j] = new double [ncols_temp[j]] {0};
		CoinCopyN(clbd_core_[j], ncols_[j], clbd_core_temp[j]);
		if (ncols_temp[j] > ncols_[j])
			CoinFill(clbd_core_temp[j] + ncols_[j], clbd_core_temp[j] + ncols_temp[j], -COIN_DBL_MAX);

		cubd_core_temp[j] = new double [ncols_temp[j]] {0};
		CoinCopyN(cubd_core_[j], ncols_[j], cubd_core_temp[j]);
		if (ncols_temp[j] > ncols_[j])
			CoinFill(cubd_core_temp[j] + ncols_[j], cubd_core_temp[j] + ncols_temp[j], COIN_DBL_MAX);
		
		ctype_core_temp[j] = new char [ncols_temp[j]] {0};
		CoinCopyN(ctype_core_[j], ncols_[j], ctype_core_temp[j]);
		if (ncols_temp[j] > ncols_[j])
			CoinFill(ctype_core_temp[j] + ncols_[j], ctype_core_temp[j] + ncols_temp[j], 'C');

		rlbd_core_temp[j] = new double [nrows_temp[j]] {0};
		CoinCopyN(rlbd_core_[j], nrows_[j], rlbd_core_temp[j]);
		rubd_core_temp[j] = new double [nrows_temp[j]] {0};
		CoinCopyN(rubd_core_[j], nrows_[j], rubd_core_temp[j]);
		if (nrows_temp[j] > nrows_[j]) 
		{
			double * startl = rlbd_core_temp[j] + nrows_[j];
			double * startu = rubd_core_temp[j] + nrows_[j];
			int length;
			
			for (i = qc_rstart[j]; i < qc_rstart[j] + nqrows_stg[j]; i++)
			{
				length = n_core[i];
				CoinFill(startl, startl + length, 0.0);
				CoinFill(startu, startu + length, 0.0);

				startl[length] = qc_row_core_->rhs[i] - 1;
				startl[length + 1] = -qc_row_core_->rhs[i] - 1;

				startu[length] = qc_row_core_->rhs[i] - 1;
				startu[length + 1] = -qc_row_core_->rhs[i] - 1;

				startl += length + 2;
				startu += length + 2;
			}
		}
	}
#ifdef DSP_DEBUG
for (j = 0; j < nstgs_; j++)
{
	cout << "obj_core_temp[" << j << "]:" << endl;
	DSPdebug(DspMessage::printArray(ncols_temp[j], obj_core_temp[j]));
}
for (j = 0; j < nstgs_; j++)
{
	cout << "clbd_core_temp[" << j << "]:" << endl;
	DSPdebug(DspMessage::printArray(ncols_temp[j], clbd_core_temp[j]));
}
for (j = 0; j < nstgs_; j++)
{
	cout << "cubd_core_temp[" << j << "]:" << endl;
	DSPdebug(DspMessage::printArray(ncols_temp[j], cubd_core_temp[j]));
}
for (j = 0; j < nstgs_; j++)
{
	cout << "rlbd_core_temp[" << j << "]:" << endl;
	DSPdebug(DspMessage::printArray(nrows_temp[j], rlbd_core_temp[j]));
}
for (j = 0; j < nstgs_; j++)
{
	cout << "rubd_core_temp[" << j << "]:" << endl;
	DSPdebug(DspMessage::printArray(nrows_temp[j], rubd_core_temp[j]));
}
#endif
	/* update clbd_core_temp and rlbd_scen_, rubd_scen_ of added variables and constraints */
	for (i = 0; i < nstgs_; i++)
	{
		int cpos = ncols_[i];
		int rpos = nrows_[i];
		
		for (j = 0; j < nqrows_stg[i]; j++)
		{
			int row_idx = qc_rstart[i] + j;
			clbd_core_temp[i][cpos + n_core[row_idx] + 1] = 0;
			
			if (i == 1)
			{
				for (s = 0; s < nscen_; s++)
				{
					if (fabs(qc_row_scen_[s]->rhs[j] - qc_row_core_->rhs[row_idx]) > 1e-8)
					{
						rlbd_scen_[s]->insert(rpos + n_core[row_idx] + 1, qc_row_scen_[s]->rhs[j] - 1);
						rlbd_scen_[s]->insert(rpos + n_core[row_idx] + 2, -qc_row_scen_[s]->rhs[j] - 1);

						rubd_scen_[s]->insert(rpos + n_core[row_idx] + 1, qc_row_scen_[s]->rhs[j] - 1);
						rubd_scen_[s]->insert(rpos + n_core[row_idx] + 2, -qc_row_scen_[s]->rhs[j] - 1);
					}
				}
			}

			cpos += n_core[row_idx] + 2;
			rpos += n_core[row_idx] + 2;
		}
	}

	/* if additional 1st stg vars are appended, we need to adjust indices for 2nd stg vars in rows_core_, mat_scen_, qc_row_core_, qc_row_scen_, indices_core, and indices_scen */
	if (ncols_temp[0] > ncols_[0])
	{
		int diff = ncols_temp[0] - ncols_[0];

		/* for 2nd stg linear constraints */
		for (i = 0; i < nrows_[1]; i++)
		{	
			for (j = 0; j < rows_core_[rstart_[1] + i]->getNumElements(); j++)
			{
				if (rows_core_[rstart_[1] + i]->getIndices()[j] >= cstart_[1])
					rows_core_[rstart_[1] + i]->getIndices()[j] += diff;
			}
			for (s = 0; s < nscen_; s++)
			{	
				for (j = mat_scen_[s]->getNumRows()-1; j >= 0; --j)
				{
					for (k = mat_scen_[s]->getNumCols()-1; k >= 0; --k)
					{
						double element = mat_scen_[s]->getCoefficient(j,k);
						if (fabs(element) > 1e-8)
						{
							if (k >= cstart_[1])
							{
								cout << "modify ind " << k  << " to " << k + diff << endl;
								mat_scen_[s]->modifyCoefficient(j, k + diff, element);
								mat_scen_[s]->modifyCoefficient(j, k, 0.0);			
							}
						}
					}
				}
			}
		}
		/* for 2nd stg quadratic constraints */
		for (i = qc_rstart[1]; i < nqrows; i++)
		{
			for (j = 0; j < qc_row_core_->linnzcnt[i]; j++)
			{
				if (qc_row_core_->linind[i][j] >= cstart_[1])
					qc_row_core_->linind[i][j] += diff;
			}
			for (j = 0; j < qc_row_core_->quadnzcnt[i]; j++)
			{
				if (qc_row_core_->quadrow[i][j] >= cstart_[1])
					qc_row_core_->quadrow[i][j] += diff;
					qc_row_core_->quadcol[i][j] += diff;
			}
			for (j = 0; j < indices_core[i].size(); j++)
			{
				if (indices_core[i][j] >= cstart_[1])
					indices_core[i][j] += diff;
			}

			for (s = 0; s < nscen_; s++)
			{
				int row_idx = i - qc_rstart[1];
				for (j = 0; j < qc_row_scen_[s]->linnzcnt[row_idx]; j++)
				{
					if (qc_row_scen_[s]->linind[row_idx][j] >= cstart_[1])
						qc_row_scen_[s]->linind[row_idx][j] += diff;
				}
				for (j = 0; j < qc_row_scen_[s]->quadnzcnt[row_idx]; j++)
				{
					if (qc_row_scen_[s]->quadrow[row_idx][j] >= cstart_[1])
						qc_row_scen_[s]->quadrow[row_idx][j] += diff;
						qc_row_scen_[s]->quadcol[row_idx][j] += diff;
				}
				for (j = 0; j < indices_scen[s][row_idx].size(); j++)
				{
					if (indices_scen[s][row_idx][j] >= cstart_[1])
						indices_scen[s][row_idx][j] += diff;
				}
			}
		}
	}

	/* copy rows_core to rows_core_temp */
	CoinPackedVector ** rows_core_temp  = new CoinPackedVector * [nrows_core_temp];
	/* copy 1st stage rows */
	for (i = 0; i < rstart_[1]; i++)
		rows_core_temp[i] = new CoinPackedVector(*(rows_core_[i]));
	/* copy 2nd stage rows */
	for (i = 0; i < nrows_[1]; i++)
		rows_core_temp[rstart_temp[1] + i] = new CoinPackedVector(*(rows_core_[rstart_[1] + i]));

	/* add affine core rows of each quadratic rows to rows_core_temp*/
	int *cstart_of_col_to_add = new int [nstgs_];
	int *rstart_of_row_to_add = new int [nstgs_];
	for (i = 0; i < nstgs_; i++)
	{
		cstart_of_col_to_add[i] = cstart_temp[i] + ncols_[i];
		rstart_of_row_to_add[i] = rstart_temp[i] + nrows_[i];
	}

	for (i = 0; i < nqrows; i++)
	{
		if (i >= qc_rstart[1])
		{
			/* if it is a 2nd stg quadratic row, then append from the current end of the 2nd stg linear constraints */
			addAffineRowOfL(L_core[i], n_core[i], indices_core[i], rows_core_temp, cstart_of_col_to_add[1], rstart_of_row_to_add[1], qc_row_core_->linnzcnt[i], qc_row_core_->linind[i], qc_row_core_->linval[i]);

			cstart_of_col_to_add[1] += n_core[i] + 2;
			rstart_of_row_to_add[1] += n_core[i] + 2;
		}
		else 
		{
			/* if it is a 1st stg quadratic row, then append from the current end of the 1st stg linear constraints */
			addAffineRowOfL(L_core[i], n_core[i], indices_core[i], rows_core_temp, cstart_of_col_to_add[0], rstart_of_row_to_add[0], qc_row_core_->linnzcnt[i], qc_row_core_->linind[i], qc_row_core_->linval[i]);

			cstart_of_col_to_add[0] += n_core[i] + 2;
			rstart_of_row_to_add[0] += n_core[i] + 2;
		}
	}
	for (i = 0; i < nstgs_; i++)
	{
		assert(cstart_of_col_to_add[i] - cstart_temp[i] == ncols_temp[i]);
		assert(rstart_of_row_to_add[i] - rstart_temp[i] == nrows_temp[i]);
	}
	FREE_ARRAY_PTR(cstart_of_col_to_add);
	FREE_ARRAY_PTR(rstart_of_row_to_add);
	
	/* This is for updating mat_scen_: 
	 * (1) add affine rows of each scenario quadratic rows to rows_scen_temp
	 * (2) compare the affine rows to construct mat_scen_ */
	for (s = 0; s < nscen_; s++)
	{	
		int cstart_of_col_to_add = cstart_temp[1] + ncols_[1];
		int rstart_of_row_to_add = 0;

		CoinPackedVector ** rows_scen_temp = new CoinPackedVector * [nrows_temp[1]];
		
		for (i = qc_rstart[1]; i < qc_rstart[1] + nqrows_stg[1]; i++)
		{
			int row_idx = i - qc_rstart[1];
			
			addAffineRowOfL(L_scen[s][row_idx], n_scen[s][row_idx], indices_scen[s][row_idx], rows_scen_temp, cstart_of_col_to_add, rstart_of_row_to_add, qc_row_scen_[s]->linnzcnt[row_idx], qc_row_scen_[s]->linind[row_idx], qc_row_scen_[s]->linval[row_idx]);

			cstart_of_col_to_add += n_scen[s][row_idx] + 2;
			rstart_of_row_to_add += n_scen[s][row_idx] + 2;
		}
		
		assert(cstart_of_col_to_add - cstart_temp[1] == ncols_temp[1]);
		assert(rstart_of_row_to_add + nrows_[1] == nrows_temp[1]);

#ifdef DSP_DEBUG
/* Compare rows_core_temp[rstart_temp[1] + nrows_[1] + j] vs rows_scen_temp[j] */
for (j = 0; j < rstart_of_row_to_add; j++)
{
	int core_ridx = rstart_temp[1] + nrows_[1] + j;
	for (k = 0; k < rows_core_temp[core_ridx]->getNumElements(); k++)
	{
		if (k == 0) {
			cout << rlbd_core_temp[1][nrows_[1] + j] << " <= ";
		}
			
		cout << rows_core_temp[core_ridx]->getElements()[k] << " " << rows_core_temp[core_ridx]->getIndices()[k] << " + ";

		if (k == rows_core_temp[core_ridx]->getNumElements() - 1) {
			cout << " <= " << rubd_core_temp[1][nrows_[1] + j] << endl;
		}
	}

	cout << j << "th add. scen " << s << " row: " << endl;
	int scen_ridx = j;
	for (k = 0; k < rows_scen_temp[scen_ridx]->getNumElements(); k++)
	{
		if (k == 0) {
			cout << "? <= ";
		}
			
		cout << rows_scen_temp[scen_ridx]->getElements()[k] << " " << rows_scen_temp[scen_ridx]->getIndices()[k] << " + ";

		if (k == rows_scen_temp[scen_ridx]->getNumElements() - 1) {
			cout << " <= ?" << endl;
		}
	}
}
#endif
		
		int naddterms = 0;
		vector<int> nchgterms(nrows_temp[1] - nrows_[1], 0); /* number of changed terms in each added linear rows */
		vector<int> inds;									 /* the set of indices for changed terms */
		vector<int> vals;									 /* the set of values for changed terms */
		
		int *row_starts = new int [nrows_temp[1]];
		assert(nrows_[1] == mat_scen_[s]->getMajorDim());
		/* copy original row_starts */
		CoinCopyN(mat_scen_[s]->getVectorStarts(), mat_scen_[s]->getMajorDim(), row_starts);
		/* calculate the last element of the original row_starts */
		row_starts[nrows_[1]] = mat_scen_[s]->getVectorStarts()[nrows_[1]-1] + mat_scen_[s]->getVectorLengths()[nrows_[1]-1];
		
		for (j = 0; j < nrows_temp[1] - nrows_[1]; j++)
		{
			int core_ridx = rstart_temp[1] + nrows_[1] + j;
			
			int * idx_core = rows_core_temp[core_ridx]->getIndices();
			double * val_core = rows_core_temp[core_ridx]->getElements();
			int nzcnt_core = rows_core_temp[core_ridx]->getNumElements();
			// double rlbd_core = rlbd_core_temp[1][nrows_[1] + j];
			// double rubd_core = rubd_core_temp[1][nrows_[1] + j];

			int scen_ridx = j;
			
			int * idx_scen = rows_scen_temp[scen_ridx]->getIndices();
			double * val_scen = rows_scen_temp[scen_ridx]->getElements();
			int nzcnt_scen = rows_scen_temp[scen_ridx]->getNumElements();
			// double rlbd_core = qc_row_scen_[s]->rhs[1];
			// double rubd_core = rubd_core_temp[1][nrows_[1] + j];

			for (k = 0; k < nzcnt_core; k++)
			{
				int idx;
				int found = false;
				for (int l = 0; l < nzcnt_scen; l++)
				{
					if (idx_scen[l] == idx_core[k])
					{
						idx = l;
						found = true;
						break;
					}
				}
				if (found)
				{
					if (fabs(val_core[k] - val_scen[idx]) > 1e-8) 
					{
#ifdef DSP_DEBUG
cout << "idx_scen: " << idx_scen[idx] << ", idx_core[k]: " << idx_core[k] << endl; 
cout << "diff: " << fabs(val_core[k] - val_scen[idx]) << endl;
#endif
						nchgterms[j]++;
						inds.push_back(idx_scen[idx]);
						vals.push_back(val_scen[idx]);
					}
				}
				else 
				{
#ifdef DSP_DEBUG					
cout << "not found: idx_core[k]: " << idx_core[k] << endl; 
#endif
					nchgterms[j]++;
					inds.push_back(idx_core[k]);
					vals.push_back(0);
				}
			}
			
			/* update row_starts */
			if (j < nrows_temp[1] - nrows_[1] - 1)
				row_starts[nrows_[1] + j + 1] = row_starts[nrows_[1] + j] + nchgterms[j];

			naddterms += nchgterms[j];
		}
		
		assert(inds.size() == naddterms);
		assert(vals.size() == naddterms);

		int num_matrix_elements = mat_scen_[s]->getNumElements() + naddterms;
		double *row_elements = new double [num_matrix_elements];
		/* copy original matrix elements to row_elements */
		CoinCopyN(mat_scen_[s]->getElements(), mat_scen_[s]->getNumElements(), row_elements);
		int *row_indices = new int [num_matrix_elements];
		/* copy original matrix indices to row_indices */
		CoinCopyN(mat_scen_[s]->getIndices(), mat_scen_[s]->getNumElements(), row_indices);
		/* append new matrix indices and elements */
		int pos = mat_scen_[s]->getNumElements();
		for (int l = 0; l < naddterms; l++)
		{
			row_elements[pos] = vals[l];
			row_indices[pos] = inds[l];
			pos++;
		}
		assert(pos == num_matrix_elements);
		/* construct lens */
		vector<int> lens(nrows_temp[1]);
		for (k = 0; k < nrows_[1]; k++)
		{
			lens[k] = mat_scen_[s]->getVectorLengths()[k];
		}
		for (k = nrows_[1]; k < nrows_temp[1]; ++k)
			lens[k] = nchgterms[k - nrows_[1]];

		CoinPackedMatrix * mat_scen_temp = new CoinPackedMatrix(false, ncols_core_temp, nrows_temp[1], num_matrix_elements, 
																	row_elements, row_indices, row_starts, &lens[0]);

#ifdef DSP_DEBUG
cout << "mat_scen_temp[" << s << "]:" << endl;
for (int l = 0; l < mat_scen_temp->getNumRows(); l++){
	cout << l << "th row: " << endl;
	if (l != mat_scen_temp->getNumRows() - 1){
		for (int k = mat_scen_temp->getVectorStarts()[l]; k < mat_scen_temp->getVectorStarts()[l+1]; k++)
			cout << mat_scen_temp->getElements()[k] << " " << mat_scen_temp->getIndices()[k] << " + ";
		cout << endl;
	} else {
		for (int k = mat_scen_temp->getVectorStarts()[l]; k < mat_scen_temp->getNumElements(); k++)
			cout << mat_scen_temp->getElements()[k] << " " << mat_scen_temp->getIndices()[k] << " + ";	
		cout << endl;
	}
	cout << endl;
}
#endif
		delete mat_scen_[s];
		mat_scen_[s] = mat_scen_temp;

		FREE_ARRAY_PTR(row_starts);
		FREE_ARRAY_PTR(row_indices);
		FREE_ARRAY_PTR(row_elements);
		FREE_ARRAY_PTR(rows_scen_temp);
	}

	/* add corresponding socp */
	for (i = 0; i < nstgs_; i++)
	{	
		/* first stage */
		vector<char> sense(nqrows_stg[i], 'L');
		vector<double> rhs(nqrows_stg[i], 0.0);
		vector<vector<int>> linind(nqrows_stg[i]);
		vector<vector<double>> linval(nqrows_stg[i]);
		vector<vector<int>> quadrow(nqrows_stg[i]);
		vector<vector<int>> quadcol(nqrows_stg[i]);
		vector<vector<double>> quadval(nqrows_stg[i]);

		int pos = cstart_temp[i] + ncols_[i];
		for (j = 0; j < nqrows_stg[i]; j++)
		{
			/* j the row */
			quadrow[j].push_back(pos + n_core[j + qc_rstart[i]] + 1);
			quadcol[j].push_back(pos + n_core[j + qc_rstart[i]] + 1);
			quadval[j].push_back(-1);

			for (k = 0; k < n_core[j + qc_rstart[i]] + 1; k++)
			{
				quadrow[j].push_back(pos + k);
				quadcol[j].push_back(pos + k);
				quadval[j].push_back(1);
			}

			pos += n_core[j + qc_rstart[i]] + 2;
		}
		
		if (i == 0)
		{
			delete qc_row_core_;
			qc_row_core_ = new QuadRowData(nqrows_stg[i], qc_rstart[i], sense, rhs, linind, linval, quadrow, quadcol, quadval);
		}
		else
		{
			for (s = 0; s < nscen_; s++)
			{ 
				delete qc_row_scen_[s];
				qc_row_scen_[s] = new QuadRowData(nqrows_stg[i], 0, sense, rhs, linind, linval, quadrow, quadcol, quadval);
			}
		}
	}

	FREE_ARRAY_PTR(nrows_);
	FREE_ARRAY_PTR(ncols_);
	FREE_ARRAY_PTR(rstart_);
	FREE_ARRAY_PTR(cstart_);
	FREE_2D_ARRAY_PTR(nstgs_, obj_core_);
	FREE_2D_ARRAY_PTR(nstgs_, clbd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, cubd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, ctype_core_);
	FREE_2D_ARRAY_PTR(nstgs_, rlbd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, rubd_core_);
	FREE_ARRAY_PTR(rows_core_);

	nrows_ = nrows_temp;
	ncols_ = ncols_temp;
	rstart_ = rstart_temp;
	cstart_ = cstart_temp;
	obj_core_ = obj_core_temp;
	clbd_core_ = clbd_core_temp;
	cubd_core_ = cubd_core_temp;
	ctype_core_ = ctype_core_temp;
	rlbd_core_ = rlbd_core_temp;
	rubd_core_ = rubd_core_temp;
	rows_core_ = rows_core_temp;

	nrows_core_ = nrows_core_temp;
	ncols_core_ = ncols_core_temp;

	FREE_ARRAY_PTR(n_core);
	FREE_2D_ARRAY_PTR(nscen_, n_scen);
	FREE_ARRAY_PTR(nqrows_stg);
	FREE_2D_ARRAY_PTR(nqrows, L_core);
	for (s = 0; s < nscen_; s++) 
		FREE_2D_ARRAY_PTR(qc_row_scen_[s]->nqrows, L_scen[s]);
	FREE_ARRAY_PTR(L_scen);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** read quadratic data file */
DSP_RTN_CODE StoModel::readQuad(const char * smps, const char * filename, bool chg_to_socp)
{
	BGN_TRY_CATCH

	map<string, int> map_varName_index;

	if (!mapVarnameIndex(map_varName_index, smps)) 
	{
		char msg[128];
		sprintf(msg, "Unable to map variables to their indices\n");
		throw msg;
	}

	char quad_cor[128], quad_sto[128], quad_tim[128];
	sprintf(quad_cor, "%s.cor", filename); 
	sprintf(quad_sto, "%s.sto", filename); 
	sprintf(quad_tim, "%s.tim", filename); 
	ifstream corfile(quad_cor);
	ifstream stofile(quad_sto);
	ifstream timfile(quad_tim);
	
	/* read core file */
	string fileName;
	string rhsName = "NULL";
	int nqrows = 0;
	map<string, int> map_qrowName_index;
	vector<char> sense;
	vector<double> rhs;
	vector<vector<int>> linind;
	vector<vector<double>> linval;
	vector<vector<int>> quadrow;
	vector<vector<int>> quadcol;
	vector<vector<double>> quadval; 

	string line;
	vector<string> items;
	map<string, int>::iterator it;

	int i, j, k;

	/* convert to 'L' sense constraint */
	if (corfile.is_open()) {

		while (getline(corfile, line, '\n')) {

			split(line, items);
			if (items[0].find("NAME") != string::npos) 
			{
				fileName = items[1];
			}
			else if (items[0].find("ROWS") != string::npos)
			{
				/** read row data */
				while(getline(corfile, line, '\n')) 
				{
					split(line, items);
					
					if (items[0].find("COLUMNS") != string::npos || items[0].find("RHS") != string::npos || items[0].find("ENDATA") != string::npos) {
						long pos = corfile.tellg();
						int leng = line.length()+1;
						corfile.seekg(pos-leng);
						break;
					}
					else if (items[0] == "G")
						sense.push_back('G');
					else if (items[0] == "L")
						sense.push_back('L');
					else {
						char msg[128];
						sprintf(msg, "Quadratic constraints sense must be 'G' or 'L'\n");
						throw msg;
					}
					map_qrowName_index[items[1]] = nqrows;
					nqrows++;
				}

				/** allocate memory */
				assert(nqrows == sense.size());
				assert(nqrows == map_qrowName_index.size());
			}
			else if (items[0].find("COLUMNS") != string::npos)
			{
				/** read linind, linval */
				linind.resize(nqrows);
				linval.resize(nqrows);
				
				while (getline(corfile, line, '\n')) 
				{	
					split(line, items);
					
					if (items[0].find("RHS") != string::npos || items[0].find("QCMATRIX") != string::npos || items[0].find("ENDATA") != string::npos) {
						long pos = corfile.tellg();
						int leng = line.length()+1;
						corfile.seekg(pos-leng);
						break;
					}

					int nterms = items.size() - 1;
					if (nterms <= 0 && nterms % 2 != 0) {
						char msg[128];
						sprintf(msg, "Each row in the COLUMNS section needs to start with a variable name and is followed by a series of pairs of a row name and its associated coefficient\n");
						throw msg;
					}

					int i = 1;
					while (i < nterms) {
						
						it = map_qrowName_index.find(items[i++]);
						if (it == map_qrowName_index.end()) 
						{
							char msg[128];
							sprintf(msg, "All rows should be declared before COLUMNS data\n");
							throw msg; 
						} 
						else 
						{
							linind[it->second].push_back(map_varName_index[items[0]]);
							if (sense[it->second] == 'G')
								linval[it->second].push_back(-stod(items[i++]));
							else 
								linval[it->second].push_back(stod(items[i++]));
						}
					}
				}
			}
			else if (items[0].find("RHS") != string::npos)
			{
				/** read rhs_ */
				rhs.resize(nqrows);
				while (getline(corfile, line, '\n')) {
					
					split(line, items);
					
					if (items[0].find("QCMATRIX") != string::npos || items[0].find("ENDATA") != string::npos || items[0].find("COLUMNS") != string::npos || items[0].find("ROWS") != string::npos){
						long pos = corfile.tellg();
						int leng = line.length()+1;
						corfile.seekg(pos-leng);
						break;
					}
					
					if (rhsName == "NULL")
						rhsName = items[0];

					it = map_qrowName_index.find(items[1]);
					if (it == map_qrowName_index.end()) 
					{
						char msg[128];
						sprintf(msg, "All rows should be declared before RHS data\n");
						throw msg; 
					} 
					else 
					{
						// rhs[it->second] = stod(items[2]);
						if (sense[it->second] == 'G')
							rhs[it->second] = -stod(items[2]);
						else 
							rhs[it->second] = stod(items[2]);
					}
				}
			}
			else if (items[0].find("QCMATRIX") != string::npos) 
			{	
				/** start reading quad term data */		
				quadrow.resize(nqrows);
				quadcol.resize(nqrows);
				quadval.resize(nqrows);	

				it = map_qrowName_index.find(items[1]);
				if (it == map_qrowName_index.end()) 
				{
					char msg[128];
					sprintf(msg, "All rows should be declared before QCMATRIX data\n");
					throw msg; 
				} 
				else
				{
					/** read quadrow_, quadcol_, quadval_ */
					while (getline(corfile, line, '\n')) {
						split(line, items);
						
						if (items[0].find("QCMATRIX") != string::npos || items[0].find("ENDATA") != string::npos || items[0].find("RHS") != string::npos){
							long pos = corfile.tellg();
							int leng = line.length()+1;
							corfile.seekg(pos-leng);
							break;
						}

						quadrow[it->second].push_back(map_varName_index[items[0]]);
						quadcol[it->second].push_back(map_varName_index[items[1]]);
						// quadval[it->second].push_back(stod(items[2]));
						if (sense[it->second] == 'G')
							quadval[it->second].push_back(-stod(items[2]));
						else 
							quadval[it->second].push_back(stod(items[2]));
					}
				}
			}
			else if (items[0].find("ENDATA") != string::npos) 
			{
				break;
			}
			else {
				char msg[128];
				sprintf(msg, "Quadratic data file encountered unexpected input: %s\n", line.c_str());
				throw msg;
			}
		}
	    corfile.close();
    } else {
        cout << "Unable to open quad core file";
        return 1;
    }

	for (i = 0; i < nqrows; i++)
		sense[i] = 'L';

	/* read time file */
	map<string, int> map_stgName_index;
	vector<int> rstart(nstgs_);	
	if (timfile.is_open()) {
		while(getline(timfile, line)) 
		{
			split(line, items);
			if (items[0].find("PERIODS") != string::npos) 
				break;
		}
		
		int nstg = 0;
		
		while(getline(timfile, line)) {
			
			if (line.find("ENDATA") != string::npos)
				break;
			
			split(line, items);
			
			map_stgName_index[items[1]] = nstg;
			it = map_qrowName_index.find(items[0]);
			if (it == map_qrowName_index.end()) 
			{
				char msg[128];
				sprintf(msg, "There is an undeclared row in time file.\n");
				throw msg; 
			} 
			else 
			{
				rstart[nstg] = it->second;
			}

			nstg++;
		}
		if (nstg != nstgs_){
			char msg[128];
			sprintf(msg, "number of stages of smps data and that of quad data do not match in time file.\n");
			throw msg; 
		}
			
		timfile.close();
    } else {
        cout << "Unable to open quad time file";
        return 1;
	}

	/* store data in StoModel */
	int nqrows_core = rstart[1] - rstart[0];
	int nqrows_scen = nqrows - rstart[1];
	if (chg_to_socp)
		qc_row_core_ = new QuadRowData(nqrows, rstart[0], sense, rhs, linind, linval, quadrow, quadcol, quadval);
	else
	{
		qc_row_core_ = new QuadRowData(nqrows_core, rstart[0], sense, rhs, linind, linval, quadrow, quadcol, quadval);
		/* check whether there is a coupling quadratic row */
		for (i = 0; i < qc_row_core_->nqrows; i++) {
			for (j = 0; j < qc_row_core_->linnzcnt[i]; j++) {
				if (qc_row_core_->linind[i][j] >= ncols_[0]) {
					char msg[128];
					sprintf(msg, "There is a second stage var in a linear term of a first stage quadratic row and chg_to_scop is turned off. If there is a coupling quadratic constraint, please turn on chg_to_scop.\n");
					throw msg; 
				}
			}
			for (j = 0; j < qc_row_core_->quadnzcnt[i]; j++) {
				if (qc_row_core_->quadcol[i][j] >= ncols_[0] || qc_row_core_->quadrow[i][j] >= ncols_[0]) {
					char msg[128];
					sprintf(msg, "There is a second stage var in a quadratic term of a first stage quadratic row and chg_to_scop is turned off. If there is a coupling quadratic constraint, please turn on chg_to_scop.\n");
					throw msg; 
				}
			}
		}
	}
	qc_row_scen_ = new QuadRowData * [nscen_];
	for (int s = 0; s < nscen_; s++) 
	{	
		qc_row_scen_[s] = new QuadRowData(nqrows_scen, rstart[1], sense, rhs, linind, linval, quadrow, quadcol, quadval);
		if (!chg_to_socp)
		{
			/* check whether there is a coupling quadratic row */
			for (i = 0; i < qc_row_scen_[s]->nqrows; i++) {
				for (j = 0; j < qc_row_scen_[s]->linnzcnt[i]; j++) {
					if (qc_row_scen_[s]->linind[i][j] < ncols_[0]) {
						char msg[128];
						sprintf(msg, "There is a first stage var in a linear term of a second stage quadratic row and chg_to_scop is turned off. If there is a coupling quadratic constraint, please turn on chg_to_scop.\n");
						throw msg; 
					}
				}
				for (j = 0; j < qc_row_scen_[s]->quadnzcnt[i]; j++) {
					if (qc_row_scen_[s]->quadcol[i][j] < ncols_[0] || qc_row_scen_[s]->quadrow[i][j] < ncols_[0]) {
						char msg[128];
						sprintf(msg, "There is a first stage var in a quadratic term of a second stage quadratic row and chg_to_scop is turned off. If there is a coupling quadratic constraint, please turn on chg_to_scop.\n");
						throw msg; 
					}
				}
			}
		}
	}	

	if (stofile.is_open()) {
		while(getline(stofile, line)) 
		{
			split(line, items);
			
			if (items[0].find("SCENARIOS") != string::npos) 
				break;
		}
		
		int nscen = 0;
		while(getline(stofile, line)) {
			
			if (line.find("ENDATA") != string::npos)
				break;
			
			split(line, items);
			
			if (items[0] == "SC") {

				it = map_stgName_index.find(items[4]);
				if (it == map_stgName_index.end()) 
				{
					char msg[128];
					sprintf(msg, "There is an undeclared stage name in stofile.\n");
					throw msg; 
				} 
				int stg = it->second;
				QuadRowData *qc = qc_row_scen_[nscen];
				nscen++;
				
				while(getline(stofile, line)) {
					split(line, items);
					
					if (items[0].find("SC") != string::npos || items[0].find("ENDATA") != string::npos){
						long pos = stofile.tellg();
						int leng = line.length()+1;
						stofile.seekg(pos-leng);
						break;
					}

					string rname;
					if (items.size() == 3) {
						/* linterm or rhs */
						it = map_qrowName_index.find(items[1]);
					} else if (items.size() == 4) {
						/* quadterm */
						it = map_qrowName_index.find(items[2]);
					} else {
						char msg[128];
						sprintf(msg, "There is an error in stofile.\n");
						throw msg; 
					}
					
					if (it == map_stgName_index.end()) 
					{
						char msg[128];
						sprintf(msg, "There is an undeclared stage name in stofile.\n");
						throw msg; 
					}
					int row_index = it->second - rstart[stg];
					bool found = false;
					if (items.size() == 3) 
					{
						if (items[0].find(rhsName) == string::npos) {
							it = map_varName_index.find(items[0]);
							if (it == map_varName_index.end()) 
							{
								char msg[128];
								sprintf(msg, "There is an undeclared variable name in stofile.\n");
								throw msg; 
							}
							/* linterm */
							for (i = 0; i < qc->linnzcnt[row_index]; i++) {
								if (qc->linind[row_index][i] == it->second) {
									qc->linval[row_index][i] = stod(items[2]);
									found = true;
									break;
								}
							}
						} else {
							qc->rhs[row_index] = stod(items[2]);
							found = true;
						}
					}
					else if (items.size() == 4) {
						/* quadterm */
						vector<int> var_index(2);
						for (j = 0; j < 2; j++) {
							it = map_varName_index.find(items[j]);
							if (it == map_varName_index.end()) 
							{
								char msg[128];
								sprintf(msg, "There is an undeclared variable name in stofile.\n");
								throw msg; 
							}
							else var_index[j] = it->second;
						}
						for (i = 0; i < qc->quadnzcnt[row_index]; i++) 
						{	
							if (qc->quadrow[row_index][i] == var_index[0] && qc->quadcol[row_index][i] == var_index[1]) 
							{
								qc->quadval[row_index][i] = stod(items[3]);
								found = true;
								break;
							}
						}
					}
					if (!found) {
						char msg[128];
						sprintf(msg, "constraint structure in stofile must agree with corfile.\n");
						throw msg; 
					}
				}
			}
		}
		if (nscen != nscen_){
			char msg[128];
			sprintf(msg, "number of scenarios of smps data and that of quad data do not match in stofile.\n");
			throw msg; 
		}
		stofile.close();
    } else {
        cout << "Unable to open quad stochastic file";
        return 1;
	}

#ifdef DSP_DEBUG
printQuadRows(-1);
for (i = 0; i < nscen_; i++)
	printQuadRows(i);	
#endif

	if (chg_to_socp)
		chgToSocp(rstart);

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

	normalizeProbability();

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
			 */
			for (int ss = 0; ss < nscen_; ++ss)
			{
				for (int j = 0; j < ncols_[1]; ++j)
				{
					wass_dist_[r][ss] += pow(fabs((*obj_scen_[s])[j] - (*obj_scen_[ss])[j]), lp_norm);
					if ((*clbd_scen_[s])[j] > -1.e+20 && (*clbd_scen_[ss])[j] > -1.e+20)
						wass_dist_[r][ss] += pow(fabs((*clbd_scen_[s])[j] - (*clbd_scen_[ss])[j]), lp_norm);
					if ((*cubd_scen_[s])[j] < 1.e+20 && (*cubd_scen_[ss])[j] < 1.e+20)
						wass_dist_[r][ss] += pow(fabs((*cubd_scen_[s])[j] - (*cubd_scen_[ss])[j]), lp_norm);
				}
				for (int i = 0; i < nrows_[1]; ++i)
				{
					if ((*rlbd_scen_[s])[i] > -1.e+20 && (*rlbd_scen_[ss])[i] > -1.e+20)
						wass_dist_[r][ss] += pow(fabs((*rlbd_scen_[s])[i] - (*rlbd_scen_[ss])[i]), lp_norm);
					if ((*rubd_scen_[s])[i] < 1.e+20 && (*rubd_scen_[ss])[i] < 1.e+20)
						wass_dist_[r][ss] += pow(fabs((*rubd_scen_[s])[i] - (*rubd_scen_[ss])[i]), lp_norm);
					for (int j = 0; j < mat_scen_[s]->getNumCols(); ++j)
					{
						wass_dist_[r][ss] += pow(fabs(mat_scen_[s]->getCoefficient(i, j) - mat_scen_[ss]->getCoefficient(i, j)), lp_norm);
					}
				}
				/* Quadratic constraints */
				if (hasQuadraticRowScenario()) {
					QuadRowData * qc_s = qc_row_scen_[s];
					QuadRowData * qc_ss = qc_row_scen_[ss];
					for (int i = 0; i < qc_s->nqrows; i++) 
					{
						wass_dist_[r][ss] += pow(fabs(qc_s->rhs[i] - qc_ss->rhs[i]), 2);

						for (int j = 0; j < qc_s->linnzcnt[i]; j++) {
							wass_dist_[r][ss] += pow(fabs(qc_s->linval[i][j] - qc_ss->linval[i][j]), 2);
						}
						for (int j = 0; j < qc_s->quadnzcnt[i]; j++) {
							wass_dist_[r][ss] += pow(fabs(qc_s->quadval[i][j] - qc_ss->quadval[i][j]), 2);
						}
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

void StoModel::normalizeProbability()
{
	double prob_sum = 0.0;
	for (int s = 0; s < nscen_; ++s)
		prob_sum += prob_[s];

	if (fabs(prob_sum - 1.0) > 1.e-8)
		for (int s = 0; s < nscen_; ++s)
			prob_[s] /= prob_sum;
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

void StoModel::copyCoreQuadraticObjective(CoinPackedMatrix *&qobj_coupling, CoinPackedMatrix *&qobj_ncoupling, int stg)
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

	if (adjustProbability && prob_[scen] > 1e-8)
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

	QuadRowData *qc;
	if (s < 0) {
	// 	qc = qc_row_core_;
	// 	cout << "Core quad constr: ";
	} else {
		qc = qc_row_scen_[s];
		cout << "Scen " << s << " quad constr: ";
	}
	for (int i = 0; i < qc->nqrows; i++) 
	{
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
