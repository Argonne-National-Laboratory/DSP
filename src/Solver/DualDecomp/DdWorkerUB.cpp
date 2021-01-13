/*
 * DdWorkerUB.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG
#include "Model/DecTssModel.h"
#include "Solver/DualDecomp/DdWorkerUB.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"

#ifdef DSP_HAS_SCIP
#include "Solver/DualDecomp/SCIPconshdlrBendersDd.h"
#endif

DdWorkerUB::DdWorkerUB(
	DecModel *model, /**< model pointer */
	DspParams *par,	 /**< parameter pointer */
	DspMessage *message /**< message pointer */) : DdWorker(model, par, message),
												   bestub_(COIN_DBL_MAX),
												   mat_mp_(NULL),
												   rlbd_org_(NULL),
												   rubd_org_(NULL),
												   osi_(NULL),
												   ub_(0.0)
{
}

DdWorkerUB::DdWorkerUB(const DdWorkerUB &rhs) : DdWorker(rhs),
												bestub_(rhs.bestub_),
												primsols_(rhs.primsols_),
												ub_(rhs.ub_)
{
	// number of subproblems to take care of
	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	// allocate memory
	mat_mp_ = new CoinPackedMatrix *[nsubprobs];
	rlbd_org_ = new double *[nsubprobs];
	rubd_org_ = new double *[nsubprobs];
	osi_ = new DspOsi *[nsubprobs];
	for (int s = 0; s < nsubprobs; ++s)
	{
		mat_mp_[s] = new CoinPackedMatrix(*(rhs.mat_mp_[s]));
		rlbd_org_[s] = new double[osi_[s]->si_->getNumRows()];
		rubd_org_[s] = new double[osi_[s]->si_->getNumRows()];
		osi_[s] = rhs.osi_[s]->clone();
		CoinCopyN(rhs.rlbd_org_[s], osi_[s]->si_->getNumRows(), rlbd_org_[s]);
		CoinCopyN(rhs.rubd_org_[s], osi_[s]->si_->getNumRows(), rubd_org_[s]);
	}
}

DdWorkerUB::~DdWorkerUB()
{
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), mat_mp_);
	FREE_2D_ARRAY_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rlbd_org_);
	FREE_2D_ARRAY_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), rubd_org_);
	FREE_2D_PTR(par_->getIntPtrParamSize("ARR_PROC_IDX"), osi_);
}

DSP_RTN_CODE DdWorkerUB::init()
{
	BGN_TRY_CATCH
	/** status */
	
	status_ = DSP_STAT_MW_CONTINUE;
	/** create problem */
	DSP_RTN_CHECK_THROW(createProblem());
	
	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DdWorkerUB::createProblem()
{
#define FREE_MEMORY            \
	FREE_PTR(mat_reco)         \
	FREE_ARRAY_PTR(clbd_reco)  \
	FREE_ARRAY_PTR(cubd_reco)  \
	FREE_ARRAY_PTR(ctype_reco) \
	FREE_ARRAY_PTR(obj_reco)

	/** recourse problem data */
	CoinPackedMatrix * mat_reco = NULL;
	double * clbd_reco   = NULL;
	double * cubd_reco   = NULL;
	double * obj_reco    = NULL;
	//CoinPackedMatrix * qobj_reco_coupling = NULL;
	//CoinPackedMatrix * qobj_reco_ncoupling = NULL;
	char *   ctype_reco  = NULL;

	/** DRO UB problem */
	double * clbd_dro = NULL;
	double * cubd_dro = NULL;
	double * obj_dro = NULL;
	double * rlbd_dro = NULL;
	double * rubd_dro = NULL;
	double * elem_dro = NULL;
	int * ind_dro = NULL;
	int * bgn_dro = NULL;
	int * len_dro = NULL;

	TssModel* tss = NULL;

	bool has_integer = false;

	BGN_TRY_CATCH

	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	try
	{
		tss = dynamic_cast<TssModel *>(model_);
	}
	catch (const std::bad_cast &e)
	{
		printf("This is not a stochastic programming problem.\n");
		return DSP_RTN_ERR;
	}

	/** allocate memory */
	mat_mp_    = new CoinPackedMatrix * [nsubprobs];
	rlbd_org_  = new double * [nsubprobs];
	rubd_org_  = new double * [nsubprobs];

	qobj_reco_coupling_ = new CoinPackedMatrix * [nsubprobs];
	qobj_reco_ncoupling_ = new CoinPackedMatrix * [nsubprobs];
	osi_       = new DspOsi * [nsubprobs];
	primsols_.resize(nsubprobs);

	for (int s = 0; s < nsubprobs; ++s)
	{

		/** copy recourse problem */
		
		if (tss->hasQuadraticObjCore() == false && tss->hasQuadraticObjScenario()==false){
			DSP_RTN_CHECK_THROW(model_->copyRecoProb(par_->getIntPtrParam("ARR_PROC_IDX")[s],
				mat_mp_[s], mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, rlbd_org_[s], rubd_org_[s]));
				qobj_reco_coupling_[s]=NULL;
				qobj_reco_ncoupling_[s]=NULL;
		}
		else{
				
			DSP_RTN_CHECK_THROW(model_->copyRecoProb(par_->getIntPtrParam("ARR_PROC_IDX")[s],
				mat_mp_[s], mat_reco, clbd_reco, cubd_reco, ctype_reco,
				obj_reco, qobj_reco_coupling_[s], qobj_reco_ncoupling_[s], rlbd_org_[s], rubd_org_[s]));
				
		}
		//printf("number of elements in qobj_reco_coupling = %d\n", qobj_reco_coupling_[s]->getNumElements());
		//PRINT_ARRAY_MSG(qobj_reco_ncoupling->getNumElements(), qobj_reco_ncoupling->getElements(), "elements in qobj_reco_ncoupling");

		for (int j = 0; j < mat_reco->getNumCols(); ++j)
			if (ctype_reco[j] != 'C')
			{
				has_integer = true;
				break;
			}

		/** creating solver interface */
		osi_[s] = createDspOsi();
		if (!osi_[s])
			throw CoinError("Failed to create DspOsi", "createProblem", "DdWorkerUB");

		/** no display */
		osi_[s]->setLogLevel(0);
		DSPdebug(osi_[s]->setLogLevel(5));
		
	    /** load problem */
	    osi_[s]->si_->loadProblem(*mat_reco, clbd_reco, cubd_reco, obj_reco, rlbd_org_[s], rubd_org_[s]);
		if (qobj_reco_ncoupling_[s]!=NULL){
			osi_[s]->loadQuadraticObjective(*qobj_reco_ncoupling_[s]);
		}
		
		for (int j = 0; j < mat_reco->getNumCols(); ++j) {
			if (ctype_reco[j] != 'C')
				osi_[s]->si_->setInteger(j);
		}
		DSPdebug(mat_reco->verifyMtx(4));
		DSPdebugMessage("number of integers: %d\n", osi_[s]->si_->getNumIntegers());

		/** allocate array size for each scenario primal solution */
		primsols_[s].resize(osi_[s]->si_->getNumCols());

		/** add quadratic constraints */
		if (tss->hasQuadraticRowCore()) 
		{
			/* nothing to do since current version only accepts noncoupling quadratic rows
			 * any coupling solution obtained from each subproblem will satisfy the first stage quadratic rows 
			 */
		}
		if (tss->hasQuadraticRowScenario()) 
		{
			QuadRowData * qc_row_data = tss->getQuaraticsRowScenario(s);

#ifdef DSP_DEBUG
			/* print qc_row_data to test whether it is successfully received or not */
			tss->printQuadRows(s);
			tss->printQuadRows(qc_row_data);
#endif
			/* Note that current version supports quadratic constraints with only second stage variables */
			int nqrow = qc_row_data->nqrows;
			int ** linind = new int * [nqrow];
			int ** quadrow = new int * [nqrow];
			int ** quadcol = new int * [nqrow];
			
			for (int i = 0; i < nqrow; i++)
			{	
				int linnzcnt = qc_row_data->linnzcnt[i];
				int quadnzcnt = qc_row_data->quadnzcnt[i];
				
				linind[i] = new int [linnzcnt];
				quadrow[i] = new int [quadnzcnt];
				quadcol[i] = new int [quadnzcnt];

				for (int j = 0; j < linnzcnt; j++) 
				{
					linind[i][j] = qc_row_data->linind[i][j] - tss->getNumCols(0);

					if (linind[i][j] < 0)
					{
						throw "Current version supports quadratic constraints in non-coupling rows only./n";
					}
				}
        
				for (int j = 0; j < quadnzcnt; j++) 
				{
					quadrow[i][j] = qc_row_data->quadrow[i][j] - tss->getNumCols(0);
					quadcol[i][j] = qc_row_data->quadcol[i][j] - tss->getNumCols(0);

					if (quadrow[i][j] < 0 || quadcol[i][j] < 0)
					{
						throw "current version support quadratic constraints in non-coupling rows only./n";
					}
				}
			}
			osi_[s]->addQuadraticRows(qc_row_data->nqrows, qc_row_data->linnzcnt, qc_row_data->quadnzcnt, qc_row_data->rhs, qc_row_data->sense, linind, qc_row_data->linval, quadrow, quadcol, qc_row_data->quadval);
			
			FREE_2D_ARRAY_PTR(nqrow, linind);
			FREE_2D_ARRAY_PTR(nqrow, quadrow);
			FREE_2D_ARRAY_PTR(nqrow, quadcol);
		}
#ifdef DSP_DEBUG
/* write in lp file to see whether the quadratic rows are successfully added to the model or not */
char filename[128];
sprintf(filename, "DdWorkerUB_scen%d", s); 
		osi_[s]->writeProb(filename, "lp");
#endif

	}

	END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

double DdWorkerUB::evaluate(int n, double *solution)
{
	std::vector<int> indices;
	std::vector<double> elements;
	for (int i = 0; i < n; ++i)
		if (fabs(solution[i]) > 1.0e-10)
		{
			indices.push_back(i);
			elements.push_back(solution[i]);
		}

	CoinPackedVector *s = new CoinPackedVector(indices.size(), &indices[0], &elements[0]);
	double ub = evaluate(s);

	delete s;

	return ub;
}

double DdWorkerUB::evaluate(CoinPackedVector *solution)
{
#define FREE_MEMORY    \
	FREE_ARRAY_PTR(Tx) \
	FREE_ARRAY_PTR(x) \
	FREE_ARRAY_PTR(xq) \
	FREE_ARRAY_PTR(xy) \
	FREE_ARRAY_PTR(qxy) \
	FREE_ARRAY_PTR(xyq) 

	double * Tx = NULL;
	double* x = NULL;
	double * xy = NULL;
	double * qxy = NULL;
	double * xyq = NULL;
	double * xq = NULL;
	TssModel* tss = NULL;

	BGN_TRY_CATCH

	tss = dynamic_cast<TssModel *>(model_);
	if (tss == NULL)
		throw "This is not a stochastic programming problem.";

	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	/** first-stage objective value */
	double cx = solution->dotProduct(tss->getObjCore(0));
	double cx_weighted = 0.0;

	/** allocate memory */
	x = solution->denseVector(tss->getNumCols(0));

	//PRINT_ARRAY_MSG(tss->getNumCols(0), x, "first stage solution");

	/** if qcqp, add quadratic part */
	if (tss->hasQuadraticObjCore()){
		if (tss->getQuadraticObjCore(0)!=NULL){
			xq = new double [tss->getNumCols(0)];
			//PRINT_ARRAY_MSG(tss->getQuadraticObjCore(0)->getNumElements(), tss->getQuadraticObjCore(0)->getElements(), "elements of x^2");
			tss->getQuadraticObjCore(0)->transposeTimes(x, xq);
			double xqx=0;
			for (int i=0; i<tss->getNumCols(0);i++){
				xqx+=x[i]*xq[i];
			}
			cx+=xqx;
			}
	}
	
	for (int s = nsubprobs - 1; s >= 0; --s) { //why from the last subprob?
		/** calculate Tx */
		Tx = new double[mat_mp_[s]->getNumRows()];
		mat_mp_[s]->times(x, Tx);

		/** adjust row bounds */
		const double *rlbd = osi_[s]->si_->getRowLower();
		const double *rubd = osi_[s]->si_->getRowUpper();
		for (int i = osi_[s]->si_->getNumRows() - 1; i >= 0; --i)
		{
			if (rlbd_org_[s][i] > -COIN_DBL_MAX)
				osi_[s]->si_->setRowLower(i, rlbd[i] - Tx[i]);
			if (rubd_org_[s][i] < COIN_DBL_MAX)
				osi_[s]->si_->setRowUpper(i, rubd[i] - Tx[i]);
		}

		assert(tss->getNumCols(1)==osi_[s]->si_->getNumCols());
		DSPdebugMessage("number of variables in the subproblem = %d\n", osi_[s]->si_->getNumCols());

		/** if QCQP, change second stage objective (coupling) if existing*/
		//if (qobj_reco_coupling_[s]->getNumElements()!=0){
			if (tss->hasQuadraticObjScenario()){
				if (qobj_reco_coupling_[s]!=NULL){
					DSPdebugMessage("there are coupling constraints\n");
				xy=new double [tss->getNumCols(0)+tss->getNumCols(1)];
				for (int i=0; i<tss->getNumCols(0)+tss->getNumCols(1);i++){
					if (i<tss->getNumCols(0)){
						xy[i]=x[i];
					}
					else{
						xy[i]=0;
					}
				}
			
				xyq = new double [tss->getNumCols(0)+tss->getNumCols(1)];
				qxy = new double [tss->getNumCols(0)+tss->getNumCols(1)];
				//PRINT_ARRAY_MSG(qobj_reco_coupling_[s]->getNumElements(), qobj_reco_coupling_[s]->getElements(), "elements in qobj_coupling");
				qobj_reco_coupling_[s]->transposeTimes(xy, xyq);
				qobj_reco_coupling_[s]->times(xy, qxy);
				vector <double> xqy;
				vector <double> ycoef;
				//PRINT_ARRAY_MSG(osi_[s]->si_->getNumCols(), osi_[s]->si_->getObjCoefficients(), "original coef of y = %f");
				for (int i=0; i< tss->getNumCols(1); i++){
					xqy.push_back(xyq[tss->getNumCols(0)+i]+qxy[tss->getNumCols(0)+i]);
					DSPdebugMessage("xqy = %f\n", xqy[i]);
					ycoef.push_back(xqy[i]+osi_[s]->si_->getObjCoefficients()[i]);
					DSPdebugMessage("coefficient of y = %f\n", ycoef[i]);
				}
				osi_[s]->si_->setObjective(&ycoef[0]);
				}

			}
		//}

		/** first-stage objective value */
		cx_weighted += cx * tss->getProbability()[par_->getIntPtrParam("ARR_PROC_IDX")[s]];

		FREE_ARRAY_PTR(Tx)

#ifdef DSP_DEBUG
      /* write in lp file to see whether the quadratic rows are successfully added to the model or not */
			char filename[128];
			sprintf(filename, "DdWorkerUB_scen_updated%d", s); 
			osi_[s]->writeProb(filename, "lp");
#endif
	}
	DSPdebugMessage("cx_weighted %e\n", cx_weighted);

	/** set initial status */
	status_ = DSP_STAT_MW_CONTINUE;

	/** solve upper bounding problem */
	DSP_RTN_CHECK_RTN_CODE(solve());

	/** update upper bound */
	if (status_ == DSP_STAT_MW_CONTINUE)
		ub_ += cx_weighted;
	else
		ub_ = COIN_DBL_MAX;

	/** restore row bounds */
	for (int s = nsubprobs - 1; s >= 0; --s)
	{
		for (int i = osi_[s]->si_->getNumRows() - 1; i >= 0; --i)
		{
			osi_[s]->si_->setRowLower(i, rlbd_org_[s][i]);
			osi_[s]->si_->setRowUpper(i, rubd_org_[s][i]);
		}
	}

	END_TRY_CATCH_RTN(FREE_MEMORY, COIN_DBL_MAX)
	FREE_MEMORY

	return ub_;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdWorkerUB::solve()
{
	double cputime;
	double walltime;

	BGN_TRY_CATCH

	double primobj = 0.0;
	double dualobj = 0.0;
	double total_cputime = 0.0;
	double total_walltime = 0.0;
	int nsubprobs = par_->getIntPtrParamSize("ARR_PROC_IDX");

	for (unsigned s = 0; s < nsubprobs; ++s)
	{
		cputime = CoinCpuTime();
		walltime = CoinGetTimeOfDay();

		/** set time limit */
		osi_[s]->setTimeLimit(
			CoinMin(CoinMax(0.01, time_remains_),
					par_->getDblParam("DD/SUB/TIME_LIM")));

		/** solve */
		osi_[s]->solve();

		/** check status. there might be unexpected results. */
		int status = osi_[s]->status();
		DSPdebugMessage("status = %d\n", status);
		switch (status)
		{
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
		case DSP_STAT_FEASIBLE:
			break;
		default:
			status_ = DSP_STAT_MW_STOP;
			message_->print(10,
							"Warning: subproblem %d solution status is %d\n", s,
							status);
			break;
		}
		if (status_ == DSP_STAT_MW_STOP)
		{
			DSPdebugMessage("status_ (dsp) = DSP_STAT_MW_STOP_, status (osi) =%d\n", status);
			primobj = COIN_DBL_MAX;
			dualobj = -COIN_DBL_MAX;
			break;
		}

		primobj += osi_[s]->getPrimObjValue();
		dualobj += osi_[s]->getDualObjValue();
		CoinCopyN(osi_[s]->si_->getColSolution(), osi_[s]->si_->getNumCols(), &primsols_[s][0]);
		total_cputime += CoinCpuTime() - cputime;
		total_walltime += CoinGetTimeOfDay() - walltime;

		/** consume time */
		time_remains_ -= CoinGetTimeOfDay() - walltime;
	}

	/** get primal objective */
	ub_ = primobj;
	DSPdebugMessage("ub_ = %e\n", ub_);
	DSPdebugMessage("status_ %d\n", status_);

	/** update statistics */
	s_statuses_.push_back(status_);
	s_primobjs_.push_back(primobj);
	s_dualobjs_.push_back(dualobj);
	s_cputimes_.push_back(total_cputime);
	s_walltimes_.push_back(total_walltime);

	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DspOsi *DdWorkerUB::createDspOsi()
{
	DspOsi *osi = NULL;
	BGN_TRY_CATCH

	switch (par_->getIntParam("DD/SUB/SOLVER")) {
	case OsiCpx:
#ifdef DSP_HAS_CPX
		osi = new DspOsiCpx();
		CPXsetintparam(dynamic_cast<DspOsiCpx *>(osi)->cpx_->getEnvironmentPtr(), CPX_PARAM_SCRIND, CPX_OFF);
#else
		throw CoinError("Cplex is not available.", "createDspOsi", "DdWorkerUB");
#endif
		break;
	case OsiGrb:
#ifdef DSP_HAS_GRB
		osi = new DspOsiGrb();
#else
		throw CoinError("Gurobi is not available.", "createDspOsi", "DdWorkerUB");
#endif
		break;
	case OsiScip:
#ifdef DSP_HAS_SCIP
		osi = new DspOsiScip();
#else
		throw CoinError("Scip is not available.", "createDspOsi", "DdWorkerUB");
#endif
		break;
	default:
		throw CoinError("Invalid parameter value", "createDspOsi", "DdWorkerUB");
		break;
	}

	END_TRY_CATCH(;)
	return osi;
}
