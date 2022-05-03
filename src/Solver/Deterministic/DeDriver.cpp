/*
 * DeDriver.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: kibaekkim
 */
// #define DSP_DEBUG

#include "DspConfig.h"
#include "Model/DecTssModel.h"
#include "Solver/Deterministic/DeDriver.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiScip.h"

DeDriver::DeDriver(
	DecModel *model, /**< model pointer */
	DspParams *par,	 /**< parameters */
	DspMessage *message /**< message pointer */) : DecSolver(model, par, message)
{
}

DeDriver::DeDriver(const DeDriver &rhs) : DecSolver(rhs) {}

DeDriver::~DeDriver()
{
}

/** initilize */
DSP_RTN_CODE DeDriver::init()
{
	BGN_TRY_CATCH

	show_copyright();

	if (model_->getFullModelNumCols() > 0) {
		primsol_.resize(model_->getFullModelNumCols());
	}

	END_TRY_CATCH(;)

	return DSP_RTN_OK;
}

/** run */
DSP_RTN_CODE DeDriver::run()
{
#define FREE_MEMORY             \
	FREE_PTR(mat)               \
	FREE_ARRAY_PTR(clbd)        \
	FREE_ARRAY_PTR(cubd)        \
	FREE_ARRAY_PTR(ctype)       \
	FREE_ARRAY_PTR(obj)         \
	FREE_PTR(qobj)              \
	FREE_ARRAY_PTR(rlbd)        \
	FREE_ARRAY_PTR(rubd)        \
	FREE_ARRAY_PTR(qc_row_scen) \
	FREE_ARRAY_PTR(linind)      \
	FREE_ARRAY_PTR(quadrow)     \
	FREE_ARRAY_PTR(quadcol)

	assert(model_);

	/** model info */
	CoinPackedMatrix *mat = NULL;
	double *clbd = NULL;
	double *cubd = NULL;
	double *obj = NULL;
	CoinPackedMatrix *qobj = NULL;
	char *ctype = NULL;
	double *rlbd = NULL;
	double *rubd = NULL;	
	
	/* quadratic row data */
	QuadRowData * qc_row_core = NULL; 
	int nscen = 0;
	QuadRowData ** qc_row_scen = NULL; 
	int ***linind = NULL;
	int ***quadrow = NULL;
	int ***quadcol = NULL;
	
	BGN_TRY_CATCH

	/** get DE model */
	DSP_RTN_CHECK_THROW(model_->getFullModel(mat, clbd, cubd, ctype, obj, qobj, rlbd, rubd));

	if (model_->isStochastic())
	{
		TssModel *tssModel;
		try
		{
			tssModel = dynamic_cast<TssModel *>(model_);
		}
		catch (const std::bad_cast &e)
		{
			printf("Model claims to be stochastic when it is not");
			return DSP_RTN_ERR;
		}

		nscen = tssModel->getNumScenarios();

		/** relax integrality? */
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[0])
		{
			for (int j = 0; j < tssModel->getNumCols(0); ++j)
			{
				ctype[j] = 'C';
			}
		}
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[1])
		{
			CoinFillN(ctype + tssModel->getNumCols(0), tssModel->getNumScenarios() * tssModel->getNumCols(1), 'C');
		}

		/** get quadratic rows data in core */
		if (tssModel->hasQuadraticRowCore()) {
			qc_row_core = tssModel->getQuaraticsRowCore();
		}

		/** get quadratic rows data for scenarios */
		if (tssModel->hasQuadraticRowScenario()) 
		{
			qc_row_scen = new QuadRowData * [nscen];
			linind = new int ** [nscen];
			quadrow = new int ** [nscen];
			quadcol = new int ** [nscen];
			
			for (int s = 0; s < nscen; s++)
			{	
				qc_row_scen[s] = tssModel->getQuaraticsRowScenario(s);
#ifdef DSP_DEBUG
				/* print qc_row_scen to test whether it is successfully received or not */
				tssModel->printQuadRows(s);
				tssModel->printQuadRows(qc_row_scen[s]);
#endif
				/* adjust indices for second stage variables */
				int nqrow = qc_row_scen[s]->nqrows;
				linind[s] = new int *[nqrow];
				quadrow[s] = new int *[nqrow];
				quadcol[s] = new int *[nqrow];

				for (int i = 0; i < nqrow; i++)
				{
					int linnzcnt = qc_row_scen[s]->linnzcnt[i];
					int quadnzcnt = qc_row_scen[s]->quadnzcnt[i];
					
					linind[s][i] = new int[linnzcnt];
					quadrow[s][i] = new int[quadnzcnt];
					quadcol[s][i] = new int[quadnzcnt];

					for (int j = 0; j < linnzcnt; j++)
					{
						if (qc_row_scen[s]->linind[i][j] > tssModel->getNumCols(0))
							linind[s][i][j] = tssModel->getNumCols(1) * s + qc_row_scen[s]->linind[i][j];
						else
							linind[s][i][j] = qc_row_scen[s]->linind[i][j];
					}

					for (int j = 0; j < quadnzcnt; j++)
					{
						if (qc_row_scen[s]->quadrow[i][j] > tssModel->getNumCols(0))
							quadrow[s][i][j] = tssModel->getNumCols(1) * s + qc_row_scen[s]->quadrow[i][j];
						else 
							quadrow[s][i][j] = qc_row_scen[s]->quadrow[i][j];

						if (qc_row_scen[s]->quadcol[i][j] > tssModel->getNumCols(0))
							quadcol[s][i][j] = tssModel->getNumCols(1) * s + qc_row_scen[s]->quadcol[i][j];
						else
							quadcol[s][i][j] = qc_row_scen[s]->quadcol[i][j];
					}
				}
			}
		}
	}
	else
	{
		if (par_->getBoolPtrParam("RELAX_INTEGRALITY")[0] ||
			par_->getBoolPtrParam("RELAX_INTEGRALITY")[1])
		{
			for (int j = 0; j < mat->getNumCols(); j++)
			{
				if (ctype[j] != 'C')
				{
					ctype[j] = 'C';
				}
			}
		}
	}

	/** create DspOsi */
	osi_ = createDspOsi();
	if (!osi_)
		throw CoinError("Failed to create DspOsi", "run", "DeDriver");

	/** set display */
	osi_->setLogLevel(par_->getIntParam("DE/SOLVER/LOG_LEVEL"));
	DSPdebug(osi_->setLogLevel(5));

	/** set number of cores */
	osi_->setNumCores(par_->getIntParam("NUM_CORES"));

	/** load problem */
	osi_->si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	//PRINT_ARRAY_MSG(qobj->getNumElements(), qobj->getElements(), "elements in qobj");
	if (qobj != NULL)
	{
		osi_->loadQuadraticObjective(*qobj);
	}
	/* add quadratic rows */
	if (model_->isStochastic()){
		if (qc_row_core) {
			osi_->addQuadraticRows(qc_row_core->nqrows, qc_row_core->linnzcnt, qc_row_core->quadnzcnt, qc_row_core->rhs, qc_row_core->sense, qc_row_core->linind, qc_row_core->linval, qc_row_core->quadrow, qc_row_core->quadcol, qc_row_core->quadval);
		}
		if (qc_row_scen) {
			for (int s = 0; s < nscen; s++) 
			{
				osi_->addQuadraticRows(qc_row_scen[s]->nqrows, qc_row_scen[s]->linnzcnt, qc_row_scen[s]->quadnzcnt, qc_row_scen[s]->rhs, qc_row_scen[s]->sense, linind[s], qc_row_scen[s]->linval, quadrow[s], quadcol[s], qc_row_scen[s]->quadval);
				FREE_2D_ARRAY_PTR(qc_row_scen[s]->nqrows, linind[s]);
				FREE_2D_ARRAY_PTR(qc_row_scen[s]->nqrows, quadrow[s]);
				FREE_2D_ARRAY_PTR(qc_row_scen[s]->nqrows, quadcol[s]);
			}
		}
	}

	for (int j = 0; j < mat->getNumCols(); j++)
	{
		if (ctype[j] != 'C')
			osi_->si_->setInteger(j);
	}

#ifdef DSP_DEBUG
	/* write in lp file to see whether the quadratic rows are successfully added to the model or not */
		char filename[128];
		sprintf(filename, "DeModel"); 
		osi_->writeProb(filename, "lp");
#endif

	/** set optimality gap tolerance */
	osi_->setRelMipGap(par_->getDblParam("DE/GAPTOL"));

	/** time limit */
	double time_limit = par_->getDblParam("DE/WALL_LIM");
	osi_->setTimeLimit(CoinMin(time_remains_, time_limit));

	/** set node limit */
	osi_->setNodeLimit(par_->getIntParam("NODE_LIM"));

	osi_->setMiqcpMethod(par_->getIntParam("DE/SOLVER/MIQCP_METHOD"));

	/** tic */
	cputime_ = CoinCpuTime();
	walltime_ = CoinGetTimeOfDay();
	
	DSPdebugMessage("start solving the deterministic model\n");
	/** solve */
	solve();
	DSPdebugMessage("solved the deterministic model\n");
	
	/** toc */
	cputime_ = CoinCpuTime() - cputime_;
	walltime_ = CoinGetTimeOfDay() - walltime_;

	/** get solutions */
	status_ = osi_->status();
	DSPdebugMessage("optimization status: %d\n", status_);
	if (status_ == DSP_STAT_OPTIMAL ||
		status_ == DSP_STAT_STOPPED_TIME ||
		status_ == DSP_STAT_STOPPED_NODE ||
		status_ == DSP_STAT_STOPPED_GAP ||
		status_ == DSP_STAT_LIM_ITERorTIME)
	{
		/** objective bounds */
		bestprimobj_ = osi_->getPrimObjValue();
		bestdualobj_ = osi_->getDualObjValue();
		DSPdebugMessage("bestprimobj_=%f, bestdualobj_=%f\n", bestprimobj_, bestdualobj_);

		/** solution */
		if (osi_->si_->getColSolution())
		{
			DSPdebugMessage("bestprimsol_=\n");
			DSPdebug(DspMessage::printArray(osi_->si_->getNumCols(), osi_->si_->getColSolution()));
			
			// make sure that the solution vector has enough space.
			if (primsol_.size() < osi_->si_->getNumCols())
				primsol_.resize(osi_->si_->getNumCols());
			CoinCopyN(osi_->si_->getColSolution(), osi_->si_->getNumCols(), &primsol_[0]);
			bestprimsol_ = primsol_;
		}
		
		/** statistics */
		if (osi_->isMip()) {
			numNodes_ = osi_->getNumNodes();
			DSPdebugMessage("numNodes_=%d\n", numNodes_);
		}
	}
	// osi_->si_->writeMps("dsp");

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

	return DSP_RTN_OK;

#undef FREE_MEMORY
}

DSP_RTN_CODE DeDriver::solve()
{
	osi_->solve();
	return DSP_RTN_OK;
}

DSP_RTN_CODE DeDriver::finalize()
{
	return DSP_RTN_OK;
}

void DeDriver::writeExtMps(const char *name)
{
#define FREE_MEMORY       \
	FREE_PTR(osi)         \
	FREE_PTR(mat)         \
	FREE_ARRAY_PTR(clbd)  \
	FREE_ARRAY_PTR(cubd)  \
	FREE_ARRAY_PTR(ctype) \
	FREE_ARRAY_PTR(obj)   \
	FREE_ARRAY_PTR(rlbd)  \
	FREE_ARRAY_PTR(rubd)

	assert(model_);

	/** model info */
	DspOsi *osi = NULL;
	CoinPackedMatrix *mat = NULL;
	double *clbd = NULL;
	double *cubd = NULL;
	double *obj = NULL;
	CoinPackedMatrix *qobj = NULL;
	char *ctype = NULL;
	double *rlbd = NULL;
	double *rubd = NULL;

	BGN_TRY_CATCH

	/** get DE model */
	DSP_RTN_CHECK_THROW(model_->getFullModel(mat, clbd, cubd, ctype, obj, qobj, rlbd, rubd));

	/** create DspOsi */
	osi = createDspOsi();
	if (!osi)
		throw CoinError("Failed to create DspOsi", "writeExtMps", "DeDriver");

	/** load problem */
	osi->si_->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	if (qobj != NULL)
	{
		osi->loadQuadraticObjective(*qobj);
	}

	for (int j = 0; j < mat->getNumCols(); j++)
	{
		if (ctype[j] != 'C')
			osi->si_->setInteger(j);
	}

	/** write mps */
	osi->writeMps(name);

	/** save memory */
	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY, ;)
#undef FREE_MEMORY
}

DspOsi *DeDriver::createDspOsi()
{
	DspOsi *osi = NULL;
	BGN_TRY_CATCH

	switch (par_->getIntParam("DE/SOLVER"))
	{
	case OsiCpx:
#ifdef DSP_HAS_CPX
		osi = new DspOsiCpx();
#else
		throw CoinError("Cplex is not available.", "createDspOsi", "DeDriver");
#endif
		break;
	case OsiGrb:
#ifdef DSP_HAS_GRB
		osi = new DspOsiGrb();
#else
		throw CoinError("Gurobi is not available.", "createDspOsi", "DeDriver");
#endif
		break;
	case OsiScip:
#ifdef DSP_HAS_SCIP
		osi = new DspOsiScip();
#else
		throw CoinError("Scip is not available.", "createDspOsi", "DeDriver");
#endif
		break;
	default:
		throw CoinError("Invalid paramter value", "createDspOsi", "DeDriver");
		break;
	}

	END_TRY_CATCH_RTN(;, osi)
	return osi;
}
