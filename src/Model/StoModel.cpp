/*
 * StoModel.cpp
 *
 *  Created on: Sep 22, 2014
 *      Author: kibaekkim
 */

#include "StoModel.h"
#include "CoinHelperFunctions.hpp"

StageData::StageData() :
		nrows_(0),
		ncols_(0),
		nints_(0),
		rstart_(0),
		cstart_(0),
		clbd_core_(NULL),
		cubd_core_(NULL),
		obj_core_(NULL),
		// qobj_core_(NULL),
		rlbd_core_(NULL),
		rubd_core_(NULL),
		ctype_core_(NULL),
		rows_core_(NULL)
		// qc_row_core_(NULL)
{
	/** nothing to do */
}

/** copy constructor */
StageData::StageData(const StageData & rhs) :
		nrows_(rhs.nrows_),
		ncols_(rhs.ncols_),
		nints_(rhs.nints_),
		rstart_(rhs.rstart_),
		cstart_(rhs.cstart_)
{
	clbd_core_ = new double [ncols_];
	cubd_core_ = new double [ncols_];
	obj_core_ = new double [ncols_];
	rlbd_core_ = new double [nrows_];
	rubd_core_ = new double [nrows_];
	ctype_core_ = new char [ncols_];
	rows_core_ = new CoinPackedVector * [nrows_];

	CoinCopyN(rhs.clbd_core_, ncols_, clbd_core_);
	CoinCopyN(rhs.cubd_core_, ncols_, cubd_core_);
	CoinCopyN(rhs.obj_core_, ncols_, obj_core_);
	CoinCopyN(rhs.rlbd_core_, nrows_, rlbd_core_);
	CoinCopyN(rhs.rubd_core_, nrows_, rubd_core_);
	CoinCopyN(rhs.ctype_core_, ncols_, ctype_core_);

	for (int i = nrows_ - 1; i >= 0; --i)
	{
		if (rhs.rows_core_[i])
			rows_core_[i] = new CoinPackedVector(*(rhs.rows_core_[i]));
		else
			rows_core_[i] = new CoinPackedVector;
	}
}
/** destructor */
StageData::~StageData()
{
	FREE_ARRAY_PTR(clbd_core_);
	FREE_ARRAY_PTR(cubd_core_);
	FREE_ARRAY_PTR(obj_core_);
	FREE_ARRAY_PTR(rlbd_core_);
	FREE_ARRAY_PTR(rubd_core_);
	FREE_ARRAY_PTR(ctype_core_);
	FREE_2D_PTR(nrows_, rows_core_);
	nrows_ = 0;
	ncols_ = 0;
	nints_ = 0;
	rstart_ = 0;
	cstart_ = 0;
}

DspScnNode::DspScnNode() :
		stg_(0),
        prob_(0.0),
		parent_(NULL),
		children_(),
		mat_scen_(NULL),
		clbd_scen_(NULL),
		cubd_scen_(NULL),
		obj_scen_(NULL),
		// qobj_scen_(NULL),
		rlbd_scen_(NULL),
		rubd_scen_(NULL)
		// qc_row_scen_(NULL),
{
	/** nothing to do */
}

/** copy constructor */
DspScnNode::DspScnNode(const DspScnNode & rhs) :
		stg_(rhs.stg_),
        prob_(rhs.prob_),
		parent_(NULL),
		children_()
{
	int ncols = sizeof(rhs.clbd_scen_) / sizeof(double);
	assert(sizeof(rhs.cubd_scen_) / sizeof(double) == ncols 
		&& sizeof(rhs.obj_scen_) / sizeof(double) == ncols);

	int nrows = sizeof(rhs.rlbd_scen_) / sizeof(double);
	assert(sizeof(rhs.rubd_scen_) / sizeof(double) == nrows);

	if (rhs.mat_scen_)
		mat_scen_ = new CoinPackedMatrix(*(rhs.mat_scen_));
	else
		mat_scen_ = new CoinPackedMatrix;

	clbd_scen_ = new double [ncols];
	cubd_scen_ = new double [ncols];
	obj_scen_ = new double [ncols];
	rlbd_scen_ = new double [nrows];
	rubd_scen_ = new double [nrows];

	CoinCopyN(rhs.clbd_scen_, ncols, clbd_scen_);
	CoinCopyN(rhs.cubd_scen_, ncols, cubd_scen_);
	CoinCopyN(rhs.obj_scen_, ncols, obj_scen_);
	CoinCopyN(rhs.rlbd_scen_, nrows, rlbd_scen_);
	CoinCopyN(rhs.rubd_scen_, nrows, rubd_scen_);
}

/** destructor */
DspScnNode::~DspScnNode()
{
	FREE_PTR(mat_scen_);
	FREE_ARRAY_PTR(clbd_scen_);
	FREE_ARRAY_PTR(cubd_scen_);
	FREE_ARRAY_PTR(obj_scen_);
	FREE_ARRAY_PTR(rlbd_scen_);
	FREE_ARRAY_PTR(rubd_scen_);
	stg_ = 0;
	prob_ = 0.0;
}

StoModel::StoModel() :
		nstgs_(0),
		nnodes_(0),
		nscen_(0),
		nrows_core_(0),
		ncols_core_(0),
		nints_core_(0),
		stage_data_(NULL),
		node_data_(),
		fromSMPS_(false)
{
	/** nothing to do */
}

/** copy constructor */
StoModel::StoModel(const StoModel & rhs) :
		nstgs_(rhs.nstgs_),
		nnodes_(rhs.nnodes_),
		nscen_(rhs.nscen_),
		nrows_core_(rhs.nrows_core_),
		ncols_core_(rhs.ncols_core_),
		nints_core_(rhs.nints_core_),
		fromSMPS_(rhs.fromSMPS_)
{
	/** allocate memory */
	stage_data_ = new StageData * [nstgs_];
	node_data_.reserve(nnodes_);

	/** copy */
	for (int s = nstgs_ - 1; s >= 0; --s)
	{
		StageData* sd = new StageData(*(rhs.stage_data_[s]));
		stage_data_[s] = sd;
	}
	
	DspScnNode* dsproot = new DspScnNode(*(rhs.node_data_[0]));
	copyAllSubNodes(*(rhs.node_data_[0]), dsproot);

	/** copy initial solutions */
	for (unsigned i = 0; i < rhs.init_solutions_.size(); ++i)
		init_solutions_.push_back(new CoinPackedVector(rhs.init_solutions_[i]));

}

void StoModel::copyAllSubNodes(const DspScnNode & rhs, DspScnNode* dsproot)
{
	node_data_.push_back(dsproot);
	for (DspScnNode * child : rhs.children_)
	{
		DspScnNode* dspchild = new DspScnNode(*(child));
		dsproot->addChild(dspchild);
		dspchild->addParent(dsproot);
		copyAllSubNodes(*(child), dspchild);
	}
}

/** destructor */
StoModel::~StoModel()
{
	FREE_2D_ARRAY_PTR(nstgs_, stage_data_);
	for (unsigned i = 0; i < node_data_.size(); ++i)
		FREE_PTR(node_data_[i]);
	for (unsigned i = 0; i < init_solutions_.size(); ++i)
		FREE_PTR(init_solutions_[i]);
	nstgs_ = 0;
	nnodes_ = 0;
	nscen_ = 0;
	nrows_core_ = 0;
	ncols_core_ = 0;
	nints_core_ = 0;
}

DSP_RTN_CODE StoModel::readSmps(const char * filename)
{	
	SmiScnModel smi;
	SmiCoreData * core = NULL;

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

    /** number of nodes */
    nnodes_ = smi.getSmiTree()->wholeTree().size();

	/** number of scenarios */
	nscen_ = smi.getNumScenarios();

	/** allocate memory */
	assert(getNumScenarios() > 0);
	stage_data_ = new StageData * [nstgs_];
    node_data_.reserve(nnodes_);

	// for (i=0; i<nstgs_;i++){
	// 	qobj_core_[i]=NULL;
	// }

	// for (i=0; i<nscen_; i++){
	// 	qobj_scen_[i]=NULL;
	// }
	/** stage information */
	nrows_core_ = core->getNumRows();
	ncols_core_ = core->getNumCols();
	nints_core_ = core->getBinaryLength() + core->getIntegerLength();

	counter = 0; /**< for checking core rows*/

	for (int s = 0; s < nstgs_; ++s) addStage(core, s);
	assert(counter == nrows_core_);

	SmiTreeNode<SmiScnNode *> * root = smi.getSmiTree()->getRoot();
	DspScnNode* dsproot = createNode(core, root->getDataPtr());
	addAllSubNodes(core, root, dsproot);

	/** mark */
	fromSMPS_ = true;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

void StoModel::addStage(
	SmiCoreData * core,
	int s
)
{
	assert(core->getNumRows(s) >= 0);
	assert(core->getNumCols(s) >= 0);
	assert(core->getRowStart(s) >= 0);
	assert(core->getColStart(s) >= 0);
	StageData* sd = new StageData();
	sd->nrows_ = core->getNumRows(s);
	sd->ncols_ = core->getNumCols(s);
	sd->nints_ = (int) core->getIntCols(s).size();
	sd->rstart_ = core->getRowStart(s);
	sd->cstart_ = core->getColStart(s);
	sd->clbd_core_ = new double [sd->ncols_];
	sd->cubd_core_ = new double [sd->ncols_];
	sd->obj_core_ = new double [sd->ncols_];
	sd->rlbd_core_ = new double [sd->nrows_];
	sd->rubd_core_ = new double [sd->nrows_];
	sd->ctype_core_ = new char [sd->ncols_];
	sd->rows_core_ = new CoinPackedVector * [sd->nrows_];
	core->copyColLower(sd->clbd_core_, s);
	core->copyColUpper(sd->cubd_core_, s);
	core->copyObjective(sd->obj_core_, s);
	core->copyRowLower(sd->rlbd_core_, s);
	core->copyRowUpper(sd->rubd_core_, s);

	/** set column types */
	CoinFillN(sd->ctype_core_, sd->ncols_, 'C');
	for (int j = 0; j < core->getBinaryLength(); ++j)
	{
		int binIndex = core->getBinaryIndices()[j];
		if (binIndex < sd->cstart_ ||
			binIndex >= sd->cstart_ + sd->ncols_)
			continue;
		sd->ctype_core_[binIndex - sd->cstart_] = 'B';
		sd->clbd_core_[binIndex - sd->cstart_] = 0.0;
		sd->cubd_core_[binIndex - sd->cstart_] = 1.0;
	}
	for (int j = 0; j < core->getIntegerLength(); ++j)
	{
		int intIndex = core->getIntegerIndices()[j];
		if (intIndex < sd->cstart_ ||
			intIndex >= sd->cstart_ + sd->ncols_)
			continue;
		sd->ctype_core_[intIndex - sd->cstart_] = 'I';
	}
	/** construct core matrix rows */
	SmiNodeData * node = core->getNode(s);
	for (int i = sd->rstart_; i < sd->rstart_ + sd->nrows_; ++i)
	{
		/**
		 * TODO: I assume the core matrix rows are well ordered.
		 */
		if (i != counter)
		{
			CoinError("Unexpected core file structure.", "readSmps", "StoModel");
		}
		
		sd->rows_core_[i-sd->rstart_] = splitRowVec(node, core, i);
		counter++;
	}

	stage_data_[s] = sd;
}

DspScnNode* StoModel::createNode(
	SmiCoreData * core,
	SmiScnNode * scnnode
)
{
	DspScnNode* nd = new DspScnNode();

	/** get stage corresponding to scenario */
	nd->stg_ = scnnode->getStage();

	/** initialize node data */
	SmiNodeData * node = scnnode->getNode();
	nd->clbd_scen_ = new double [stage_data_[nd->stg_]->ncols_];
	nd->cubd_scen_ = new double [stage_data_[nd->stg_]->ncols_];
	nd->obj_scen_ = new double [stage_data_[nd->stg_]->ncols_];
	nd->rlbd_scen_ = new double [stage_data_[nd->stg_]->nrows_];
	nd->rubd_scen_ = new double [stage_data_[nd->stg_]->nrows_];
	node->copyColLower(nd->clbd_scen_);
	node->copyColUpper(nd->cubd_scen_);
	/** update columns for binary variables since some solvers can't process properly*/
	for (int j = 0; j < stage_data_[nd->stg_]->ncols_; ++j)
	{
		if (stage_data_[nd->stg_]->ctype_core_[j] == 'B')
		{
			nd->clbd_scen_[j] = 0.0;
			nd->cubd_scen_[j] = 1.0;
		}
	}
	node->copyObjective(nd->obj_scen_);
	node->copyRowLower(nd->rlbd_scen_);
	node->copyRowUpper(nd->rubd_scen_);

	int nrows_ = stage_data_[nd->stg_]->nrows_;
	int ncols_ = stage_data_[nd->stg_]->ncols_;
	int rstart_ = stage_data_[nd->stg_]->rstart_;
	CoinPackedVectorBase ** rows_scen_ = new CoinPackedVectorBase * [nrows_];
	for (int i = rstart_; i < rstart_ + nrows_; ++i)
	{
		CoinPackedVector* node_row_copy = splitRowVec(node, core, i);
		rows_scen_[i-rstart_] = node->combineWithCoreRow(stage_data_[nd->stg_]->rows_core_[i-rstart_],node_row_copy);
		FREE_PTR(node_row_copy);
	}
	nd->mat_scen_ = new CoinPackedMatrix(false, 0.0, 0.0);
	nd->mat_scen_->setDimensions(0, ncols_);
	nd->mat_scen_->appendRows(nrows_,rows_scen_);
	FREE_2D_ARRAY_PTR(nrows_, rows_scen_);

	/** probability */
	nd->prob_ = scnnode->getProb();

	return nd;
}

void StoModel::addAllSubNodes(
	SmiCoreData * core,
	SmiTreeNode<SmiScnNode *> * root,
	DspScnNode* dsproot
)
{
	node_data_.push_back(dsproot);
	std::vector<SmiTreeNode<SmiScnNode *> *> * children = root->getChildren();
	for (SmiTreeNode<SmiScnNode *> * child : *children)
	{
		DspScnNode* dspchild = createNode(core, child->getDataPtr());
		dsproot->addChild(dspchild);
		dspchild->addParent(dsproot);
		addAllSubNodes(core, child, dspchild);
	}
	delete children;
}

const double * StoModel::getObjNode(DspScnNode* node) 
{
    return node->obj_scen_;
}

/** construct a map that maps variable names to their indices */
//fixme - instead of index, map to pair of stage and index?
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

void StoModel::setProbability(DspScnNode* node, double probability)
{
	node->prob_ = probability;
}

void StoModel::setSolution(int size, double * solution)
{
	if (size > 0)
		init_solutions_.push_back(new CoinPackedVector(size, solution));
}

// DSP_RTN_CODE StoModel::setWassersteinAmbiguitySet(double lp_norm, double eps)
// {
// 	if (nstgs_ == 0)
// 	{
// 		std::cerr << "No stochastic problem is loaded." << std::endl;
// 		return DSP_RTN_ERR;
// 	}
// 	if (ncols_[0] == 0)
// 	{
// 		std::cerr << "No stochastic programming variable is loaded." << std::endl;
// 		return DSP_RTN_ERR;
// 	}
// 	if (nscen_ == 0)
// 	{
// 		std::cerr << "No scenario is loaded." << std::endl;
// 		return DSP_RTN_ERR;
// 	}

// 	wass_eps_ = pow(eps, lp_norm);
// 	isdro_ = true;

// 	normalizeProbability();

// 	/** Count the number of reference scenarios.
// 	 * We take all scenarios with nonzero probabilities as the references (a.k.a. 
// 	 * empirical observations).
// 	 */
// 	nrefs_ = 0;
// 	for (int s = 0; s < nscen_; ++s)
// 	{
// 		if (prob_[s] > 0)
// 			nrefs_++;
// 	}

// 	refs_probability_ = new double[nrefs_];
// 	wass_dist_ = new double *[nrefs_];

// 	/** Compute the Wasserstein distances of lp_norm,
// 	 * and assign the reference probabilities.
// 	 */
// 	for (int s = 0, r = 0; s < nscen_; ++s)
// 	{
// 		if (prob_[s] > 0)
// 		{
// 			refs_probability_[r] = prob_[s];
// 			wass_dist_[r] = new double[nscen_];
// 			CoinZeroN(wass_dist_[r], nscen_);
// 			assert(mat_scen_[s]->getNumRows() == nrows_[1]);

// 			/** Compute the distances.
// 			 * 
// 			 * TODO: Can we do in parallel?
// 			 * The relevant issues need addressed first: 
// 			 * - https://github.com/kibaekkim/DSPopt.jl/issues/14
// 			 */
// 			for (int ss = 0; ss < nscen_; ++ss)
// 			{
// 				for (int j = 0; j < ncols_[1]; ++j)
// 				{
// 					wass_dist_[r][ss] += pow(fabs((*obj_scen_[s])[j] - (*obj_scen_[ss])[j]), 2);
// 					if ((*clbd_scen_[s])[j] > -1.e+20 && (*clbd_scen_[ss])[j] > -1.e+20)
// 						wass_dist_[r][ss] += pow(fabs((*clbd_scen_[s])[j] - (*clbd_scen_[ss])[j]), 2);
// 					if ((*cubd_scen_[s])[j] < 1.e+20 && (*cubd_scen_[ss])[j] < 1.e+20)
// 						wass_dist_[r][ss] += pow(fabs((*cubd_scen_[s])[j] - (*cubd_scen_[ss])[j]), 2);
// 				}
// 				for (int i = 0; i < nrows_[1]; ++i)
// 				{
// 					if ((*rlbd_scen_[s])[i] > -1.e+20 && (*rlbd_scen_[ss])[i] > -1.e+20)
// 						wass_dist_[r][ss] += pow(fabs((*rlbd_scen_[s])[i] - (*rlbd_scen_[ss])[i]), 2);
// 					if ((*rubd_scen_[s])[i] < 1.e+20 && (*rubd_scen_[ss])[i] < 1.e+20)
// 						wass_dist_[r][ss] += pow(fabs((*rubd_scen_[s])[i] - (*rubd_scen_[ss])[i]), 2);
// 					for (int j = 0; j < mat_scen_[s]->getNumCols(); ++j)
// 					{
// 						wass_dist_[r][ss] += pow(fabs(mat_scen_[s]->getCoefficient(i, j) - mat_scen_[ss]->getCoefficient(i, j)), 2);
// 					}
// 				}
// 				/* Quadratic constraints */
// 				if (hasQuadraticRowScenario()) {
// 					QuadRowData * qc_s = qc_row_scen_[s];
// 					QuadRowData * qc_ss = qc_row_scen_[ss];
// 					for (int i = 0; i < qc_s->nqrows; i++) 
// 					{
// 						wass_dist_[r][ss] += pow(fabs(qc_s->rhs[i] - qc_ss->rhs[i]), 2);

// 						for (int j = 0; j < qc_s->linnzcnt[i]; j++) {
// 							wass_dist_[r][ss] += pow(fabs(qc_s->linval[i][j] - qc_ss->linval[i][j]), 2);
// 						}
// 						for (int j = 0; j < qc_s->quadnzcnt[i]; j++) {
// 							wass_dist_[r][ss] += pow(fabs(qc_s->quadval[i][j] - qc_ss->quadval[i][j]), 2);
// 						}
// 					}
// 				}
// 				wass_dist_[r][ss] = pow(wass_dist_[r][ss], lp_norm / 2.0);
// 			}
// 			r++;
// 		}
// 	}

// 	/** scaling vector */
// 	double scaling_constant = pow(wass_eps_, 2);
// 	for (int r = 0; r < nrefs_; ++r)
// 		for (int s = 0; s < nscen_; ++s)
// 			scaling_constant += pow(wass_dist_[r][s], 2);
// 	scaling_constant = sqrt(scaling_constant);
// 	wass_eps_ /= scaling_constant;
// 	for (int r = 0; r < nrefs_; ++r)
// 		for (int s = 0; s < nscen_; ++s)
// 			wass_dist_[r][s] /= scaling_constant;

// 	/** Quadratic equations
// 	 * TODO: The quadratic objective function and constraints need to be considered.
// 	 * - https://github.com/Argonne-National-Laboratory/DSP/issues/155
// 	 */

// 	printf("[DRO] Set %d reference scenarios.\n", nrefs_);
// 	printf("[DRO] Computed the Wasserstein distances of order %f.\n", lp_norm);

// 	return DSP_RTN_OK;
// }

// void StoModel::normalizeProbability()
// {
// 	double prob_sum = 0.0;
// 	for (int s = 0; s < nscen_; ++s)
// 		prob_sum += prob_[s];

// 	if (fabs(prob_sum - 1.0) > 1.e-8)
// 		for (int s = 0; s < nscen_; ++s)
// 			prob_[s] /= prob_sum;
// }

CoinPackedVector * StoModel::splitRowVec(
		SmiNodeData * node,
		SmiCoreData * core,
		int i /**< row index */
)
{
	int* newrowids = new int [node->getRowLength(i)];
	const int* oldrowids = node->getRowIndices(i);
	for (int j = 0; j < node->getRowLength(i); ++j) 
	{
		newrowids[j] = oldrowids[j] - core->getColStart(node->getStage());
	}
	CoinPackedVector * split = new CoinPackedVector(
			node->getRowLength(i),
			newrowids,
			node->getRowElements(i));
	FREE_ARRAY_PTR(newrowids);
	return split;
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
	split->reserve(stage_data_[stg]->ncols_);

	/** insert elements */
	for (j = 0; j < stage_data_[stg]->rows_core_[i]->getNumElements(); ++j)
	{
		ind = stage_data_[stg]->rows_core_[i]->getIndices()[j] - stage_data_[stg]->cstart_;
		if (ind >= 0 && ind < stage_data_[stg]->ncols_)
		{
			split->insert(ind, stage_data_[stg]->rows_core_[i]->getElements()[j]);
		}
	}

	return split;
}

/** copy core column lower bounds */
void StoModel::copyCoreColLower(double * clbd, int stg)
{
	CoinCopyN(stage_data_[stg]->clbd_core_, stage_data_[stg]->ncols_, clbd);
}

/** copy core column upper bounds */
void StoModel::copyCoreColUpper(double * cubd, int stg)
{
	CoinCopyN(stage_data_[stg]->cubd_core_, stage_data_[stg]->ncols_, cubd);
}

/** copy core objective coefficients */
void StoModel::copyCoreObjective(double * obj, int stg)
{
	for (int j = 0; j < stage_data_[stg]->ncols_; ++j)
	{
		obj[j] = stage_data_[stg]->obj_core_[j];
	}
}

// void StoModel::copyCoreQuadraticObjective(CoinPackedMatrix *&qobj_coupling, CoinPackedMatrix *&qobj_ncoupling, int stg)
// {
// 	vector<int> colidx;
// 	vector<int> rowidx;
// 	vector<double> elements;
// 	int numq=qobj_core_[stg]->getNumElements();
// 	const CoinBigIndex * start = qobj_core_[stg]->getVectorStarts();
// 	vector<int> colidx_coupling;
// 	vector<int> rowidx_coupling;
// 	vector<double> elements_coupling;
// 	int numq_coupling = 0;
// 	vector<int> colidx_ncoupling;
// 	vector<int> rowidx_ncoupling;
// 	vector<double> elements_ncoupling;
// 	int numq_ncoupling = 0;

// 	int rowcount=0;
// 	for (int i=0; i<qobj_core_[stg]->getMajorDim();i++){
// 		if (start[i+1]-start[i]>0){
// 				for (int k=0; k<start[i+1]-start[i]; k++){
// 					rowidx.push_back(rowcount);
// 				}
// 		}
// 		rowcount++;
// 	}

// 	for (int i=0; i<numq; i++){
// 		colidx.push_back(qobj_core_[stg]->getIndices()[i]);
// 		elements.push_back(qobj_core_[stg]->getElements()[i]);
// 	}

// 	// security check
// 	assert(numq==rowidx.size());
// 	assert(numq==colidx.size());
// 	assert(numq==elements.size());
// 	for (int i=0; i<numq; i++){
// 		if (rowidx[i]>=ncols_[stg-1] && colidx[i]>=ncols_[stg-1]){
// 			rowidx_ncoupling.push_back(rowidx[i]-ncols_[stg-1]);
// 			colidx_ncoupling.push_back(colidx[i]-ncols_[stg-1]);
// 			elements_ncoupling.push_back(elements[i]);
// 			numq_ncoupling++;
// 		}
// 		else{
// 			rowidx_coupling.push_back(rowidx[i]);
// 			colidx_coupling.push_back(colidx[i]);
// 			elements_coupling.push_back(elements[i]);
// 			numq_coupling++;
// 		}
// 	}

// 	if(numq_ncoupling!=0){
// 		qobj_ncoupling = new CoinPackedMatrix(false, &rowidx_ncoupling[0], &colidx_ncoupling[0], &elements_ncoupling[0], numq_ncoupling);
// 	}
// 	else{
// 		qobj_ncoupling = NULL;
// 	}
	
// 	if(numq_coupling!=0){
// 		qobj_coupling = new CoinPackedMatrix(false, &rowidx_coupling[0], &colidx_coupling[0], &elements_coupling[0], numq_coupling);
// 	}
// 	else{
// 		qobj_coupling = NULL;
// 	}
// }

/** copy core column types */
void StoModel::copyCoreColType(char * ctype, int stg)
{
	CoinCopyN(stage_data_[stg]->ctype_core_, stage_data_[stg]->ncols_, ctype);
}

/** copy core row lower bounds */
void StoModel::copyCoreRowLower(double * rlbd, int stg)
{
	CoinCopyN(stage_data_[stg]->rlbd_core_, stage_data_[stg]->nrows_, rlbd);
}

/** copy core row upper bounds */
void StoModel::copyCoreRowUpper(double * rubd, int stg)
{
	CoinCopyN(stage_data_[stg]->rubd_core_, stage_data_[stg]->nrows_, rubd);
}

/** combine random matrix row for a given scenario */
void StoModel::combineRandRowVec(
		CoinPackedVector * row, /**< core row vector */
		int i,                  /**< row index */
		DspScnNode* node        /**< node */)
{
	int j;
	int pos = 0; /** position to sparse row vector elements */

	/** random matrix row info */
	CoinBigIndex   start = node->mat_scen_->getVectorStarts()[i];
	const int *    ind   = node->mat_scen_->getIndices() + start;
	const double * elem  = node->mat_scen_->getElements() + start;

	for (j = 0; j < node->mat_scen_->getVectorSize(i); ++j)
	{
		//fixme
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
		DspScnNode* node        /**< node */)
{
	int j, col;
	int pos = 0; /** position to sparse row vector elements */

	/** random matrix row info */
	CoinBigIndex   start = node->mat_scen_->getVectorStarts()[i];
	const int *    ind   = node->mat_scen_->getIndices() + start;
	const double * elem  = node->mat_scen_->getElements() + start;

	for (j = 0; j < node->mat_scen_->getVectorSize(i); ++j)
	{
		//fixme
		col = ind[j] - stage_data_[stg]->cstart_;
		if (col >= 0 && col < stage_data_[stg]->ncols_)
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
void StoModel::combineRandColLower(double * clbd, DspScnNode* node)
{
	for (int i = 0; i < stage_data_[node->stg_]->ncols_; ++i)
	{
		clbd[i] = node->clbd_scen_[i];
	}
}

/** combine random column upper bounds */
void StoModel::combineRandColUpper(double * cubd, DspScnNode* node)
{
	for (int i = 0; i < stage_data_[node->stg_]->ncols_; ++i)
	{
		cubd[i] = node->cubd_scen_[i];
	}
}

/** combine random objective coefficients */
void StoModel::combineRandObjective(double * obj, DspScnNode* node, double adjustProbability = 1.0)
{
	int i, j;
	for (i = 0; i < stage_data_[node->stg_]->ncols_; ++i)
	{
		obj[j] = node->obj_scen_[i] * adjustProbability;
	}
}

// void StoModel::combineRandQuadraticObjective(CoinPackedMatrix * &qobj_coupling, CoinPackedMatrix * &qobj_ncoupling, int stg, int scen, bool adjustProbability)
// {
// 	int i, j, s;
// 	int numelements = qobj_scen_[scen]->getNumElements();
// 	const int * colindex=qobj_scen_[scen]->getIndices();
// 	const CoinBigIndex * start = qobj_scen_[scen]->getVectorStarts();
// 	const int dim=qobj_scen_[scen]->getMajorDim();
// 	const double * elements=qobj_scen_[scen]->getElements();

// 	int * rowindex;
// 	rowindex = new int [numelements];
	
// 	std::vector<int> newrowindex_coupling, newcolindex_coupling;
// 	std::vector<double> newelements_coupling;

// 	std::vector<int> newrowindex_ncoupling, newcolindex_ncoupling;
// 	std::vector<double> newelements_ncoupling;

	
// 	int idx=0;
// 	/** convert start into row index */
// 	for (int k=0; k<dim; k++){
// 		for (i=0; i<start[k+1]-start[k];i++){
// 			rowindex[idx]=k;
// 			idx++;
// 		}
// 	}

// 	int sidx=stg-1;
// 	if (sidx<0) sidx=0;
// 	for (i = 0; i < numelements; ++i)
// 	{
// 		if (colindex[i]>=ncols_[sidx] && rowindex[i]>=ncols_[sidx]){
// 			newrowindex_ncoupling.push_back(rowindex[i]-ncols_[sidx]);
// 			newcolindex_ncoupling.push_back(colindex[i]-ncols_[sidx]);
// 			if (adjustProbability){
// 				newelements_ncoupling.push_back(elements[i]*prob_[scen]);
// 			}
// 			else{
// 				newelements_ncoupling.push_back(elements[i]);
// 			}
// 		}
// 		else{
// 			newrowindex_coupling.push_back(rowindex[i]);
// 			newcolindex_coupling.push_back(colindex[i]);
// 			if (adjustProbability){
// 				newelements_coupling.push_back(elements[i]*prob_[scen]);
// 			}
// 			else{
// 				newelements_coupling.push_back(elements[i]);
// 			}
// 		}
// 	}
// 	//printf("newcolindex_coupling.size() = %d\n", newcolindex_coupling.size());
// 	//for (int m=0; m<newcolindex_ncoupling.size(); m++){
// 	//	printf("newcolindex_ncoupling[%d]=%d\n", m, newcolindex_ncoupling[m]);
// 	//}
// 	if (newelements_coupling.size()==0){
// 		qobj_coupling=NULL;
// 	}
// 	else{
// 		qobj_coupling=new CoinPackedMatrix(false, &newcolindex_coupling[0], &newrowindex_coupling[0], &newelements_coupling[0], newelements_coupling.size());
// 	}
// 	if (newelements_ncoupling.size()==0){
// 		qobj_coupling=NULL;
// 	}
// 	else{
// 		qobj_ncoupling=new CoinPackedMatrix(false, &newcolindex_ncoupling[0], &newrowindex_ncoupling[0], &newelements_ncoupling[0], newelements_ncoupling.size());
// 	}
// 	//qobj_ncoupling=new CoinPackedMatrix(false, &newcolindex_ncoupling[0], &newrowindex_ncoupling[0], &newelements_ncoupling[0], newelements_ncoupling.size());
// 	//PRINT_ARRAY_MSG(newcolindex_ncoupling.size(), &newcolindex_ncoupling[0], "elements in newcolindex_ncoupling");
// 	//PRINT_ARRAY_MSG(qobj_ncoupling->getNumElements(), qobj_ncoupling->getElements(), "elements in qobj_ncoupling");

// 	delete[] rowindex;
// }

/** combine random row lower bounds */
void StoModel::combineRandRowLower(double * rlbd, DspScnNode* node)
{
	for (int i = 0; i < stage_data_[node->stg_]->nrows_; ++i)
	{
		rlbd[i] = node->rlbd_scen_[i];
	}
}

/** combine random row upper bounds */
void StoModel::combineRandRowUpper(double * rubd, DspScnNode* node)
{
	for (int i = 0; i < stage_data_[node->stg_]->nrows_; ++i)
	{
		rubd[i] = node->rubd_scen_[i];
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
	PRINT_ARRAY_OBJECT_MSG(nstgs_, *stage_data_, nrows_, "nrows_")
	PRINT_ARRAY_OBJECT_MSG(nstgs_, *stage_data_, ncols_, "ncols_")
	PRINT_ARRAY_OBJECT_MSG(nstgs_, *stage_data_, rstart_, "rstart_")
	PRINT_ARRAY_OBJECT_MSG(nstgs_, *stage_data_, cstart_, "cstart_")

	printf("\n### Core data ###\n\n");
	printf("nrows_core_ %d\n", nrows_core_);
	printf("ncols_core_ %d\n", ncols_core_);
#if 1
	for (int s = 0; s < nstgs_; ++s)
	{
		for (int i = 0; i < stage_data_[s]->nrows_; ++i)
		{
			sprintf(tmpstr, "Core matrix row %d", i);
			PRINT_COIN_PACKED_VECTOR_MSG((*(stage_data_[s]->rows_core_[i])), tmpstr)
		}
	}
#endif
	for (int t = 0; t < nstgs_; ++t)
	{
		printf("\n### Stage %d core data ###\n\n", t);
		PRINT_ARRAY_MSG(stage_data_[t]->ncols_, stage_data_[t]->clbd_core_, "clbd_core_")
		PRINT_ARRAY_MSG(stage_data_[t]->ncols_, stage_data_[t]->cubd_core_, "cubd_core_")
		PRINT_ARRAY_MSG(stage_data_[t]->ncols_, stage_data_[t]->obj_core_, "obj_core_")
		PRINT_ARRAY_MSG(stage_data_[t]->ncols_, stage_data_[t]->ctype_core_, "ctype_core_")
		PRINT_ARRAY_MSG(stage_data_[t]->nrows_, stage_data_[t]->rlbd_core_, "rlbd_core_")
		PRINT_ARRAY_MSG(stage_data_[t]->nrows_, stage_data_[t]->rubd_core_, "rubd_core_")
		// printf("\n### quadratic objective coefficient ###\n\n");
		// if (qobj_core_[t])
		// {
		// 	printf("isColOrdered %d\n", qobj_core_[t]->isColOrdered());
		// 	PRINT_ARRAY_MSG(qobj_core_[t]->getMajorDim(), qobj_core_[t]->getVectorStarts(), "VectorStarts")
		// 	PRINT_SPARSE_ARRAY_MSG(
		// 			qobj_core_[t]->getNumElements(),
		// 			qobj_core_[t]->getIndices(),
		// 			qobj_core_[t]->getElements(),
		// 			"Elements")
		// }
	}

	for (int s = 0; s < nnodes_; ++s)
	{
		printf("\n### Scenario %d data ###\n\n", s);
		printf("probability %E\n", node_data_[s]->prob_);
		printf("stage map %d\n", node_data_[s]->stg_);
#if 1
		printf("=== BEGINNING of CoinPackedMatrix mat_scen_[%d] ===\n", s);
		if (node_data_[s]->mat_scen_)
		{
			printf("isColOrdered %d\n", node_data_[s]->mat_scen_->isColOrdered());
			PRINT_ARRAY_MSG(node_data_[s]->mat_scen_->getMajorDim(), node_data_[s]->mat_scen_->getVectorStarts(), "VectorStarts")
			PRINT_SPARSE_ARRAY_MSG(
					node_data_[s]->mat_scen_->getNumElements(),
					node_data_[s]->mat_scen_->getIndices(),
					node_data_[s]->mat_scen_->getElements(),
					"Elements")
		}
		printf("=== END of CoinPackedMatrix mat_scen_[%d] ===\n", s);

		// printf("=== BEGINNING of CoinPackedMatrix qobj_scen_[%d] ===\n", s);
		// if (qobj_scen_[s])
		// {
		// 	printf("isColOrdered %d\n", qobj_scen_[s]->isColOrdered());
		// 	PRINT_ARRAY_MSG(qobj_scen_[s]->getMajorDim(), qobj_scen_[s]->getVectorStarts(), "VectorStarts")
		// 	PRINT_SPARSE_ARRAY_MSG(
		// 			qobj_scen_[s]->getNumElements(),
		// 			qobj_scen_[s]->getIndices(),
		// 			qobj_scen_[s]->getElements(),
		// 			"Elements")
		// }
		// printf("=== END of CoinPackedMatrix qobj_scen_[%d] ===\n", s);
#endif
//		PRINT_COIN_PACKED_VECTOR_MSG((*clbd_scen_[s]), "clbd_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*cubd_scen_[s]), "cubd_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*obj_scen_[s]), "obj_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*rlbd_scen_[s]), "rlbd_scen_")
//		PRINT_COIN_PACKED_VECTOR_MSG((*rubd_scen_[s]), "rubd_scen_")
	}

	printf("\n### END of printing StoModel data ###\n");
}

// DSP_RTN_CODE StoModel::printQuadRows(const int s)
// {
// 	BGN_TRY_CATCH

// 	QuadRowData *qc;
// 	if (s < 0) {
// 		qc = qc_row_core_;
// 		cout << "Core quad constr: ";
// 	} else {
// 		qc = qc_row_scen_[s];
// 		cout << "Scen " << s << " quad constr: ";
// 	}
// 	for (int i = 0; i < qc->nqrows; i++) 
// 	{
// 		for (int lt = 0; lt < qc->linnzcnt[i]; lt++)
// 		{
// 			cout << qc->linval[i][lt] << " x" << qc->linind[i][lt] << " + ";
// 		}
// 		for (int qt = 0; qt < qc->quadnzcnt[i]-1; qt++)
// 		{
// 			cout << qc->quadval[i][qt] << " x" << qc->quadrow[i][qt] << " x" << qc->quadcol[i][qt] << " + ";
// 		}
// 		cout << qc->quadval[i][qc->quadnzcnt[i]-1] << " x" << qc->quadrow[i][qc->quadnzcnt[i]-1] << " x" << qc->quadcol[i][qc->quadnzcnt[i]-1];
// 		if (qc->sense[i] == 'L')
// 			cout << " <= " << qc->rhs[i] << endl;
// 		else 
// 			cout << " >= " << qc->rhs[i] << endl;
// 	}

// 	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

// 	return DSP_RTN_OK;
// }

// DSP_RTN_CODE StoModel::printQuadRows(const QuadRowData *qc)
// {
// 	BGN_TRY_CATCH

// 	for (int i = 0; i < qc->nqrows; i++) 
// 	{
// 		cout << i << "th quad constr: ";

// 		for (int lt = 0; lt < qc->linnzcnt[i]; lt++)
// 		{
// 			cout << qc->linval[i][lt] << " x" << qc->linind[i][lt] << " + ";
// 		}
// 		for (int qt = 0; qt < qc->quadnzcnt[i]-1; qt++)
// 		{
// 			cout << qc->quadval[i][qt] << " x" << qc->quadrow[i][qt] << " x" << qc->quadcol[i][qt] << " + ";
// 		}
// 		cout << qc->quadval[i][qc->quadnzcnt[i]-1] << " x" << qc->quadrow[i][qc->quadnzcnt[i]-1] << " x" << qc->quadcol[i][qc->quadnzcnt[i]-1];
// 		if (qc->sense[i] == 'L')
// 			cout << " <= " << qc->rhs[i] << endl;
// 		else 
// 			cout << " >= " << qc->rhs[i] << endl;
// 	}

// 	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

// 	return DSP_RTN_OK;
// }

// double StoModel::getWassersteinDist(int i, int j) {
// 	if (!isdro_) return 0.0;
// 	if (i >= nrefs_) {
// 		printf("Reference index (%d) is out of range (%d).\n", i, nrefs_);
// 		return 0.0;
// 	}
// 	if (j >= nscen_) {
// 		printf("Scenario index (%d) is out of range (%d).\n", j, nscen_);
// 		return 0.0;
// 	}
// 	return wass_dist_[i][j];
// }

// double StoModel::getReferenceProbability(int i) {
// 	if (refs_probability_ == NULL) return 0.0;
// 	if (!isdro_) return 0.0;
// 	if (i >= nrefs_) {
// 		printf("Reference index (%d) is out of range (%d).\n", i, nrefs_);
// 		return 0.0;
// 	}
// 	return refs_probability_[i];
// }
