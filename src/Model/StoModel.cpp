/*
 * StoModel.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: kibaekkim
 */

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
		rlbd_core_(NULL),
		rubd_core_(NULL),
		ctype_core_(NULL),
		prob_(NULL),
		mat_scen_(NULL),
		clbd_scen_(NULL),
		cubd_scen_(NULL),
		obj_scen_(NULL),
		rlbd_scen_(NULL),
		rubd_scen_(NULL),
		priorities_(NULL),
		fromSMPS_(false)
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
		fromSMPS_(rhs.fromSMPS_)
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
	rlbd_core_  = new double * [nstgs_];
	rubd_core_  = new double * [nstgs_];
	ctype_core_ = new char * [nstgs_];
	prob_       = new double [nscen_];
	mat_scen_   = new CoinPackedMatrix * [nscen_];
	clbd_scen_  = new CoinPackedVector * [nscen_];
	cubd_scen_  = new CoinPackedVector * [nscen_];
	obj_scen_   = new CoinPackedVector * [nscen_];
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
		if (rhs.rlbd_scen_[i])
			rlbd_scen_[i] = new CoinPackedVector(*(rhs.rlbd_scen_[i]));
		else
			rlbd_scen_[i] = new CoinPackedVector;
		if (rhs.rubd_scen_[i])
			rubd_scen_[i] = new CoinPackedVector(*(rhs.rubd_scen_[i]));
		else
			rubd_scen_[i] = new CoinPackedVector;
	}
	priorities_ = new int [nints_core_];
	CoinCopyN(rhs.priorities_, nints_core_, priorities_);
}

StoModel::~StoModel()
{
	FREE_ARRAY_PTR(nrows_);
	FREE_ARRAY_PTR(ncols_);
	FREE_ARRAY_PTR(nints_);
	FREE_ARRAY_PTR(rstart_);
	FREE_ARRAY_PTR(cstart_);
	FREE_2D_PTR(nstgs_, rows_core_);
	FREE_2D_ARRAY_PTR(nstgs_, clbd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, cubd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, obj_core_);
	FREE_2D_ARRAY_PTR(nstgs_, rlbd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, rubd_core_);
	FREE_2D_ARRAY_PTR(nstgs_, ctype_core_);
	FREE_ARRAY_PTR(prob_);
	FREE_2D_PTR(nscen_, mat_scen_);
	FREE_2D_PTR(nscen_, clbd_scen_);
	FREE_2D_PTR(nscen_, cubd_scen_);
	FREE_2D_PTR(nscen_, obj_scen_);
	FREE_2D_PTR(nscen_, rlbd_scen_);
	FREE_2D_PTR(nscen_, rubd_scen_);
	scen2stg_.clear();
	nscen_ = 0;
	nstgs_ = 0;
	nrows_core_ = 0;
	ncols_core_ = 0;
	FREE_ARRAY_PTR(priorities_);
//	for (unsigned int i = 0; i < branchingHyperplanes_.size(); ++i)
//		FREE_PTR(branchingHyperplanes_[i].first);
}

STO_RTN_CODE StoModel::readSmps(const char * filename)
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
	rlbd_core_  = new double * [nstgs_];
	rubd_core_  = new double * [nstgs_];
	ctype_core_ = new char * [nstgs_];
	prob_       = new double [nscen_];
	mat_scen_   = new CoinPackedMatrix * [nscen_];
	clbd_scen_  = new CoinPackedVector * [nscen_];
	cubd_scen_  = new CoinPackedVector * [nscen_];
	obj_scen_   = new CoinPackedVector * [nscen_];
	rlbd_scen_  = new CoinPackedVector * [nscen_];
	rubd_scen_  = new CoinPackedVector * [nscen_];

	/** stage information */
	nrows_core_ = 0;
	ncols_core_ = 0;
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
		nints_core_ = 0;
		CoinFillN(ctype_core_[i], ncols_[i], 'C');
		for (j = 0; j < core->getBinaryLength(); ++j)
		{
			if (core->getBinaryIndices()[j] < cstart_[i] ||
				core->getBinaryIndices()[j] >= cstart_[i] + ncols_[i])
				continue;
			ctype_core_[i][core->getBinaryIndices()[j] - cstart_[i]] = 'B';
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

	/** branching priority */
//	if (nints_core_ > 0)
//	{
//		priorities_ = new int [nints_core_];
//		CoinFillN(priorities_, nints_core_, 1000);
//	}

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

	for (i = 0; i < nscen_; ++i)
	{
		/** get stage corresponding to scenario */
		int stg = smi.getLeafNode(i)->getStage();

		/** add mapping */
		scen2stg_.insert(std::pair<int,int>(i,stg));

		/** allocate memory */
		mat_scen_[i]  = new CoinPackedMatrix(false, 0, 0);
		mat_scen_[i]->setDimensions(0, ncols_[stg]);

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
		for (j = rstart_[stg]; j < rstart_[stg] + nrows_[stg]; ++j)
		{
			mat_scen_[i]->appendRow(
					node->getRowLength(j), node->getRowIndices(j), node->getRowElements(j));
		}
	}

	/** mark */
	fromSMPS_ = true;

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** set branching priorities */
void StoModel::setPriorities(
		int * priorities /**< the size of the array should be
		                    at least the number of integer variables in core */)
{
	if (nints_core_ > 0)
	{
		if (priorities_ == NULL)
			priorities_ = new int [nints_core_];
		CoinCopyN(priorities, nints_core_, priorities_);
	}
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
	}

	for (int s = 0; s < nscen_; ++s)
	{
		printf("\n### Scenario %d data ###\n\n", s);
		printf("probability %E\n", prob_[s]);
		printf("stage map %d\n", scen2stg_[s]);
#if 1
		printf("=== BEGINNING of CoinPackedMatrix mat_scen_[%d] ===\n", s);
		printf("isColOrdered %d\n", mat_scen_[s]->isColOrdered());
		PRINT_ARRAY_MSG(mat_scen_[s]->getMajorDim(), mat_scen_[s]->getVectorStarts(), "VectorStarts")
		PRINT_SPARSE_ARRAY_MSG(
				mat_scen_[s]->getNumElements(),
				mat_scen_[s]->getIndices(),
				mat_scen_[s]->getElements(),
				"Elements")
		printf("=== END of CoinPackedMatrix mat_scen_[%d] ===\n", s);
#endif
		PRINT_COIN_PACKED_VECTOR_MSG((*clbd_scen_[s]), "clbd_scen_")
		PRINT_COIN_PACKED_VECTOR_MSG((*cubd_scen_[s]), "cubd_scen_")
		PRINT_COIN_PACKED_VECTOR_MSG((*obj_scen_[s]), "obj_scen_")
		PRINT_COIN_PACKED_VECTOR_MSG((*rlbd_scen_[s]), "rlbd_scen_")
		PRINT_COIN_PACKED_VECTOR_MSG((*rubd_scen_[s]), "rubd_scen_")
	}

	printf("\n### END of printing StoModel data ###\n");
}
