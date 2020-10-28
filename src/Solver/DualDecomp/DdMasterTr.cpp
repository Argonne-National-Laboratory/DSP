/*
 * DdMasterTr.cpp
 *
 *  Created on: Feb 12, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#include "CoinWarmStartBasis.hpp"
#include "SolverInterface/DspOsiClp.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "SolverInterface/DspOsiOoqp.h"
#include "SolverInterface/DspOsiOoqpEps.h"
#include "Solver/DualDecomp/DdMasterTr.h"
#include "Model/TssModel.h"

DdMasterTr::DdMasterTr(
		DecModel *   model,   /**< model pointer */
		DspParams *  par,     /**< parameter pointer */
		DspMessage * message /**< message pointer */):
DdMaster(model, par, message),
nthetas_(0),
nlambdas_(0),
nus_(0),
nPs_(0),
obj_reco_(NULL),
stability_param_(0.0),
stability_center_(NULL),
trcnt_(0),
numIters_(0),
cputime_elapsed_(0.0),
walltime_elapsed_(0.0),
isSolved_(false),
cuts_(NULL),
ncuts_minor_(0),
cutdel_param_(0.5),
linerr_(0.0),
parTr_(true),
parTrSize_(0.0),
parTrDecrease_(true),
parNumCutsPerIter_(1),
parMasterAlgo_(IPM_Feasible),
parLogLevel_(0),
nstalls_(0) {
    DSPdebug(message_->logLevel_ = 9999;);
}

DdMasterTr::DdMasterTr(const DdMasterTr& rhs) :
DdMaster(rhs),
nthetas_(rhs.nthetas_),
nlambdas_(rhs.nlambdas_),
nus_(rhs.nus_),
nPs_(rhs.nPs_),
obj_reco_(rhs.obj_reco_),
stability_param_(rhs.stability_param_),
trcnt_(rhs.trcnt_),
numIters_(rhs.numIters_),
cputime_elapsed_(rhs.cputime_elapsed_),
walltime_elapsed_(rhs.walltime_elapsed_),
isSolved_(rhs.isSolved_),
cuts_age_(rhs.cuts_age_),
possiblyDelete_(rhs.possiblyDelete_),
masterobjsAtCutAdd_(rhs.masterobjsAtCutAdd_),
ncuts_minor_(rhs.ncuts_minor_),
cutdel_param_(rhs.cutdel_param_),
linerr_(rhs.linerr_),
parTr_(rhs.parTr_),
parTrSize_(rhs.parTrSize_),
parTrDecrease_(rhs.parTrDecrease_),
parNumCutsPerIter_(rhs.parNumCutsPerIter_),
parMasterAlgo_(rhs.parMasterAlgo_),
parLogLevel_(rhs.parLogLevel_),
nstalls_(rhs.nstalls_) {
	if (rhs.stability_center_) {
		stability_center_ = new double [nlambdas_+nus_+nPs_];
		CoinCopyN(rhs.stability_center_, nlambdas_+nus_+nPs_, stability_center_);
	}
	cuts_ = new OsiCuts(*(rhs.cuts_));
}

DdMasterTr::~DdMasterTr()
{
	FREE_2D_ARRAY_PTR(model_->getNumSubproblems(), obj_reco_);
	FREE_ARRAY_PTR(stability_center_);
	FREE_PTR(cuts_);
}

/** initialize */
DSP_RTN_CODE DdMasterTr::init()
{
	BGN_TRY_CATCH

	DdMaster::init();

	/** read parameters */
	parTr_ = par_->getBoolParam("DD/TR");
	parTrSize_ = par_->getDblParam("DD/TR/SIZE");
	parTrDecrease_ = par_->getBoolParam("DD/TR/DECREASE");
	parNumCutsPerIter_ = par_->getIntParam("DD/NUM_CUTS_PER_ITER");
#ifndef DSP_HAS_OOQP
	if (par_->getIntParam("DD/MASTER_ALGO")  == IPM_Feasible)
		throw CoinError("Invalid parameter value", "init", "DdMasterTr");
#endif
	parMasterAlgo_ = par_->getIntParam("DD/MASTER_ALGO");
	if (model_->isDro() && parMasterAlgo_ == IPM_Feasible)
	{
		printf("-- DRO cannot use IPM_Feasible option.\n"
			   "-- The master problem will use IPM instead.\n");
		parMasterAlgo_ = IPM;
	}
	parLogLevel_ = par_->getIntParam("LOG_LEVEL");
	DSPdebugMessage("Trust region size %f\n", parTrSize_);

	/** create problem */
	DSP_RTN_CHECK_THROW(createProblem());

	/** clock */
	cputime_elapsed_  = CoinCpuTime();
	walltime_elapsed_ = CoinGetTimeOfDay();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterTr::solve()
{
	BGN_TRY_CATCH

	double cputime  = CoinCpuTime();
	double walltime = CoinGetTimeOfDay();

	/** trust-region may trigger infeasibility */
	bool resolve = true;

	while (resolve) {
		/** solve */
		getSiPtr()->resolve();
	
		/** mark as solved */
		isSolved_ = true;
	
		/** solver status */
		status_ = osi_->status();
		switch(status_)
		{
		case DSP_STAT_PRIM_INFEASIBLE:
			DSPdebugMessage("The master is infeasible and increases the trust region.\n");
			// increase the trust-region size
			stability_param_ *= 2;
			setTrustRegion(stability_param_, stability_center_);
			break;
		case DSP_STAT_OPTIMAL:
		case DSP_STAT_LIM_ITERorTIME:
		case DSP_STAT_STOPPED_GAP:
		case DSP_STAT_STOPPED_NODE:
		case DSP_STAT_STOPPED_TIME:
		{
			/** objective value */
			primobj_ = osi_->getPrimObjValue();
			/** get solution */
			CoinCopyN(getSiPtr()->getColSolution(), getSiPtr()->getNumCols(), &primsol_[0]);
#ifdef DSP_DEBUG
			printf("Master solution (obj %+e):\n", primobj_);
			DspMessage::printArray(getSiPtr()->getNumCols(), &primsol_[0]);
#endif
	
			/** update statistics */
			s_statuses_.push_back(status_);
			s_primobjs_.push_back(osi_->getPrimObjValue());
			s_dualobjs_.push_back(osi_->getDualObjValue());
			double * s_primsol = new double [getSiPtr()->getNumCols()];
			CoinCopyN(getSiPtr()->getColSolution(), getSiPtr()->getNumCols(), s_primsol);
			s_primsols_.push_back(s_primsol);
			s_primsol = NULL;
			s_cputimes_.push_back(CoinCpuTime() - cputime);
			s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);
			message_->print(3, "Master solution time %.2f sec.\n", CoinGetTimeOfDay() - walltime);

			resolve = false;
			DSPdebug(getSiPtr()->writeMps("master"));
	
			break;
		}
		default:
			message_->print(0, "Warning: master solution status is %d\n", status_);
			status_ = DSP_STAT_MW_STOP;
			resolve = false;
			DSPdebug(getSiPtr()->writeMps("master"));
			break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterTr::createProblem()
{
#define FREE_MEMORY        \
	FREE_ARRAY_PTR(ctype); \
	FREE_ARRAY_PTR(clbd);  \
	FREE_ARRAY_PTR(cubd);  \
	FREE_ARRAY_PTR(obj);   \
	FREE_ARRAY_PTR(rlbd);  \
	FREE_ARRAY_PTR(rubd);  \
	FREE_ARRAY_PTR(bgn);   \
	FREE_ARRAY_PTR(len);   \
	FREE_ARRAY_PTR(ind);   \
	FREE_ARRAY_PTR(elem);  \
	FREE_PTR(mat);

	int i, j, k, pos;
	int ncols, nrows, nzcnt;
	char * ctype = NULL;
	double * clbd = NULL;
	double * cubd = NULL;
	double * obj  = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;
	CoinBigIndex * bgn = NULL;
	int * len = NULL;
	int * ind = NULL;
	double * elem = NULL;
	CoinPackedMatrix * mat = NULL;

	TssModel *tss = NULL;

	BGN_TRY_CATCH

	if (model_->isStochastic()) {
		try {
			tss = dynamic_cast<TssModel*>(model_);
		} catch (const std::bad_cast &e) {
			printf("Error: Model claims to be stochastic when it is not\n");
			return DSP_RTN_ERR;
		}
	}

	/** Dual variables for SO */
	nthetas_  = model_->getNumSubproblems();//CoinMin(model_->getNumSubproblems(), parNumCutsPerIter_);
	nlambdas_ = model_->getNumCouplingRows();

	/** Additional dual variables for DRO */
	if (model_->isDro()) {
		nus_ = model_->getNumSubproblems() * model_->getNumReferences();
		nPs_ = model_->getNumSubproblems();
	}

	/** LP dimension */
	if (model_->nonanticipativity())
	{
		/** initial normalization constraint for nonanticipativity constraints */
		nrows = model_->getNumCouplingCols(); // number of first-stage variables
		nzcnt = nrows * model_->getNumSubproblems();

		/** Additional constraints for DRO */
		if (model_->isDro()) {
			nrows += 1 + model_->getNumReferences() + model_->getNumSubproblems();
			for (i = 0; i < model_->getNumSubproblems(); ++i) {
				for (j = 0; j < model_->getNumReferences(); ++j) {
					if (fabs(model_->getWassersteinDist(j,i)) < 1e-10) continue;
					nzcnt++;
				}
			}
			nzcnt += nus_ * 2 + nPs_;
		}
	}
	else
	{
		nrows = 0;
		nzcnt = 0;
	}
	ncols = nthetas_ + nlambdas_;

	/** Increase number of columns for DRO */
	ncols += nus_ + nPs_;
	DSPdebugMessage("nrows %d ncols %d nzcnt %d nthetas_ %d nlambdas_ %d\n",
                    nrows, ncols, nzcnt, nthetas_, nlambdas_);

	if (model_->isDro()) {
		obj_reco_ = new double * [tss->getNumScenarios()];
		for (int j = 0; j < tss->getNumScenarios(); ++j)
			tss->copyRecoObj(j, obj_reco_[j], false);
	}

	/** allocate memory */
	ctype = new char [ncols];
	clbd  = new double[ncols];
	cubd  = new double[ncols];
	obj   = new double[ncols];
	rlbd  = new double[nrows];
	rubd  = new double[nrows];
	bgn   = new CoinBigIndex[nrows + 1];
	len   = new int[nrows];
	ind   = new int[nzcnt];
	elem  = new double[nzcnt];

	/** all continuous variables */
	CoinFillN(ctype, ncols, 'C');

	/** c */
	CoinFillN(obj, nthetas_, 1.0);
	CoinZeroN(obj + nthetas_, ncols - nthetas_);

	/** trust region */
	stability_param_ = parTr_ ? parTrSize_ : COIN_DBL_MAX;
	if (parTr_)
	{
		stability_center_ = new double [ncols-nthetas_];
		if (model_->nonanticipativity()) {
			for (i = 0; i < model_->getNumCouplingCols(); ++i) {
				for (j = 0; j < model_->getNumSubproblems(); ++j) {
					k = model_->getNumCouplingCols()*j+i;
					stability_center_[k] = model_->getCouplingColsObjs()[i] / model_->getNumSubproblems();
				}
			}

			/** stability center for DRO model */
			if (model_->isDro())
			{
				for (i = 0; i < model_->getNumReferences(); ++i) {
					for (j = 0; j < model_->getNumSubproblems(); ++j) {
						if (i == j)
							stability_center_[nlambdas_+model_->getNumReferences()*j+i] = model_->getReferenceProbability(i);
						else
							stability_center_[nlambdas_+model_->getNumReferences()*j+i] = 0.0;
					}
				}
				for (i = 0; i < nPs_; ++i) {
					if (i < model_->getNumReferences())
						stability_center_[nlambdas_+nus_+i] = model_->getReferenceProbability(i);
					else
						stability_center_[nlambdas_+nus_+i] = 0.0;
				}
			}
		} else {
			CoinZeroN(stability_center_, ncols-nthetas_);
		}
	}

	/** bounds */
	CoinFillN(clbd, nthetas_, -COIN_DBL_MAX);
	CoinFillN(cubd, nthetas_, +COIN_DBL_MAX);

	/** trust region bound */
	if (model_->nonanticipativity()) {
		for (i = 0; i < model_->getNumCouplingCols(); ++i) {
			double avg_c = model_->getCouplingColsObjs()[i] / model_->getNumSubproblems();
			for (j = 0; j < model_->getNumSubproblems(); ++j) {
				k = nthetas_+model_->getNumCouplingCols()*j+i;
				clbd[k] = avg_c - stability_param_;
				cubd[k] = avg_c + stability_param_;
			}
		}

		/** create row-ordered constraint matrix */
		pos = 0;
		int rpos = 0;
		// lambdas
		for (i = 0; i < model_->getNumCouplingCols(); ++i)
		{
			rlbd[rpos] = model_->getCouplingColsObjs()[i];
			rubd[rpos] = model_->getCouplingColsObjs()[i];

			bgn[rpos] = pos;
			for (j = 0; j < model_->getNumSubproblems(); ++j)
			{
				assert(nthetas_ + j * model_->getNumCouplingCols() + i < ncols);
				ind[pos] = nthetas_ + j * model_->getNumCouplingCols() + i;
				elem[pos] = 1.0;
				pos++;
			}
			len[rpos] = pos - bgn[rpos];
			rpos++;
		}

		if (model_->isDro()) {

			/** additional column bounds for DRO */
			for (i = 0; i < model_->getNumReferences(); ++i) {
				for (j = 0; j < model_->getNumSubproblems(); ++j) {
					clbd[nthetas_+nlambdas_+model_->getNumReferences()*j+i] = CoinMax(0.0, stability_center_[nlambdas_+model_->getNumReferences()*j+i] - stability_param_);
					cubd[nthetas_+nlambdas_+model_->getNumReferences()*j+i] = CoinMin(1.0, stability_center_[nlambdas_+model_->getNumReferences()*j+i] + stability_param_);
				}
			}
			for (i = 0; i < nPs_; ++i) {
				clbd[nthetas_+nlambdas_+nus_+i] = CoinMax(0.0, stability_center_[nlambdas_+nus_+i] - stability_param_);
				cubd[nthetas_+nlambdas_+nus_+i] = CoinMin(1.0, stability_center_[nlambdas_+nus_+i] + stability_param_);
			}

			/** additional constraints for DRO */

			// Wasserstein distance
			rlbd[rpos] = -COIN_DBL_MAX;
			rubd[rpos] = model_->getWassersteinSize();

			bgn[rpos] = pos;
			for (i = 0; i < model_->getNumSubproblems(); ++i) {
				for (j = 0; j < model_->getNumReferences(); ++j) {
					if (fabs(model_->getWassersteinDist(j,i)) < 1e-10) continue;
					ind[pos] = nthetas_ + nlambdas_ + i * model_->getNumReferences() + j;
					elem[pos] = model_->getWassersteinDist(j,i);
					pos++;
				}
			}
			len[rpos] = pos - bgn[rpos];
			rpos++;

			// u's
			for (i = 0; i < model_->getNumReferences(); ++i) {
				rlbd[rpos] = model_->getReferenceProbability(i);
				rubd[rpos] = model_->getReferenceProbability(i);

				bgn[rpos] = pos;
				for (j = 0; j < model_->getNumSubproblems(); ++j) {
					ind[pos] = nthetas_ + nlambdas_ + j * model_->getNumReferences() + i;
					elem[pos] = 1.0;
					pos++;
				}
				len[rpos] = pos - bgn[rpos];
				rpos++;
			}

			// P's
			for (j = 0; j < model_->getNumSubproblems(); ++j) {
				rlbd[rpos] = 0;
				rubd[rpos] = 0;

				bgn[rpos] = pos;
				for (i = 0; i < model_->getNumReferences(); ++i) {
					ind[pos] = nthetas_ + nlambdas_ + j * model_->getNumReferences() + i;
					elem[pos] = 1.0;
					pos++;
				}
				ind[pos] = nthetas_ + nlambdas_ + + nus_ + j;
				elem[pos] = -1.0;
				pos++;
				len[rpos] = pos - bgn[rpos];
				rpos++;
			}
		}

		assert(rpos == nrows);
		assert(pos == nzcnt);
		bgn[nrows] = pos;

	} else {
		CoinFillN(clbd+nthetas_, nlambdas_, -COIN_DBL_MAX);
		CoinFillN(cubd+nthetas_, nlambdas_, +COIN_DBL_MAX);

		/** nonnegative or nonpositive multipliers according to sense */
		for (i = 0; i < nlambdas_; i++)
		{
			if (model_->getSenseCouplingRow(i) == 'L')
				clbd[nthetas_ + i] = 0;
			else if (model_->getSenseCouplingRow(i) == 'G')
				cubd[nthetas_ + i] = 0;
		}
	}

	for (int j = nthetas_; j < ncols; ++j) {
		DSPdebugMessage("j = %d, clbd = %e, stability_center_ = %e, cubd = %e\n", j, clbd[j], stability_center_[j-nthetas_], cubd[j]);
		assert(clbd[j] <= stability_center_[j-nthetas_]);
		assert(cubd[j] >= stability_center_[j-nthetas_]);
	}
#ifdef DSP_DEBUG
	DSPdebugMessage("stability_center_:\n");
	DspMessage::printArray(ncols-nthetas_, stability_center_);
#endif

	/** constraint matrix */
	mat = new CoinPackedMatrix(false, ncols, nrows, nzcnt, elem, ind, bgn, len);
	// DSPdebug(mat->verifyMtx(4));

	/** create solver interface */
	if (parMasterAlgo_ != IPM_Feasible) {
		osi_ = createDspOsi();
		if (!osi_) throw CoinError("Failed to create DspOsi", "createProblem", "DdMasterTr");
	}

	DSPdebugMessage("parMasterAlgo_ %d\n", parMasterAlgo_);
	switch (parMasterAlgo_) {
	case Simplex:
		osi_->use_simplex();
		break;
	case IPM:
		osi_->use_barrier();
		break;
	case IPM_Feasible:

#ifdef DSP_HAS_OOQP
		osi_ = new DspOsiOoqpEps();
#else
		throw CoinError("DspOsiOoqpEps is not available.", "createProblem", "DdMasterTr");
#endif
		break;
	default:
		throw CoinError("Invalid parameter value", "createProblem", "DdMasterTr");
		break;
	}

	if (!osi_) throw CoinError("Failed to create DspOsi", "createProblem", "DdMasterTr");
	DSPdebugMessage("Created master algorithm\n");

	osi_->setNumCores(par_->getIntParam("NUM_CORES"));

	/** [MAX]imization */
	getSiPtr()->setObjSense(-1);

	/** copy problem data */
	getSiPtr()->loadProblem(*mat, clbd, cubd, obj, rlbd, rubd);
	DSPdebugMessage("Loaded problem data\n");

	/** allocate memory for solution */
	primsol_.resize(ncols);
	CoinFillN(&primsol_[0], nthetas_, COIN_DBL_MAX);
	if (model_->nonanticipativity()) {
		CoinCopyN(stability_center_, ncols-nthetas_, &primsol_[nthetas_]);
	} else {
		CoinZeroN(&primsol_[nthetas_], nlambdas_);
	}

	if (model_->isStochastic()) {
		bestdualsol_.resize(ncols-nthetas_);
	}

	/** initialize cut pool */
	cuts_ = new OsiCuts;

	/** set print level */
	osi_->setLogLevel(CoinMax(0,par_->getIntParam("LOG_LEVEL")-1));

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY;

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMasterTr::updateProblem()
{
	BGN_TRY_CATCH

	int nCutsAdded = 0;
	/** current primal objective value */
	double curprimobj = 0.0;
	if (isSolved_)
		curprimobj = osi_->getPrimObjValue();

	/** calculate primal/dual objectives */
	double newprimal = 0.0;
	double newdual = 0.0;
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
	{
		DSPdebugMessage("subdualobj_[%d] = %e\n", s, subdualobj_[s]);
		newprimal += subprimobj_[s];
		//newdual += subdualobj_[s];
		newdual += subprimobj_[s];
	}

	/** update trust region FIRST, because this does not change problem. */
	if (parTr_)
	{
		DSPdebugMessage("stability_param_ %e\n", stability_param_);
		DSPdebugMessage("newdual %+e, bestdualobj %+e, curprimobj %+e\n", newdual, bestdualobj_, curprimobj);
		DSPdebugMessage("Serious test %+e >= 0\n", newdual - bestdualobj_ - 1.0e-4 * (curprimobj - bestdualobj_));
		if (newdual >= bestdualobj_ + 1.0e-4 * (curprimobj - bestdualobj_))
		{
			message_->print(2, "TR  %s STEP: dual objective %e", isSolved_ ? "SERIOUS" : "INITIAL", newdual);

			/** mark cuts not to be deleted */
			for (int i = cuts_->sizeCuts() - nCutsAdded; i < cuts_->sizeCuts(); ++i)
				possiblyDelete_[i] = false;

			if (isSolved_)
			{
				/** update proximal point */
				CoinCopyN(&primsol_[nthetas_], getSiPtr()->getNumCols()-nthetas_, stability_center_);
				message_->print(3, ", updated proximal point");

				/** possibly delete cuts */
				//DSP_RTN_CHECK_THROW(possiblyDeleteCuts(newprimal));

				/** is solution boundary? */
				if (isSolutionBoundary() &&
					newdual >= bestdualobj_ + 0.5 * (curprimobj - bestdualobj_))
				{
					/** increase trust region */
					stability_param_ = CoinMin(2. * stability_param_, 1.0e+4);
					message_->print(3, ", increased trust region size %e", stability_param_);
				}

				/** set trust region */
				setTrustRegion(stability_param_, stability_center_);
			}
			else {
				nCutsAdded = addCuts();
				ncuts_minor_ += nCutsAdded;
			}

			/** update dual bound */
			bestdualobj_ = newdual;
			trcnt_ = 0;
			nstalls_ = 0;

			/** update best solution */
			DSPdebugMessage("getSiPtr()->getNumCols()-nthetas_ [%d] == bestdualsol_.size() [%d]\n", 
				getSiPtr()->getNumCols()-nthetas_, (int) bestdualsol_.size());
			assert(getSiPtr()->getNumCols()-nthetas_==bestdualsol_.size());
			CoinCopyN(&primsol_[nthetas_], getSiPtr()->getNumCols()-nthetas_, &bestdualsol_[0]);

			message_->print(2, "\n");
		}
		else
		{
			/** add cuts and increase minor cut counter */
			nCutsAdded = addCuts();
			ncuts_minor_ += nCutsAdded;

			message_->print(4, "TR  master has %d rows and %d cols after adding %d cuts.\n",
						getSiPtr()->getNumRows() + nCutsAdded, getSiPtr()->getNumCols(), nCutsAdded);

			/** null step */
			message_->print(3, "TR  null step: dual objective %e", newdual);

			/** The following rule is from Linderoth and Wright (2003) */
			double rho = CoinMin(1.0, stability_param_) * CoinMax(bestdualobj_ - newdual, linerr_) / (curprimobj - bestdualobj_);
			if (rho > 0) trcnt_++;
			if (rho >= 3 || (trcnt_ >= 3 && fabs(rho - 2.) < 1.0))
			{
				/** decrease trust region */
				stability_param_ = CoinMax(1.0e-2, stability_param_/CoinMin(rho, 4.));
				message_->print(3, ", decreased trust region size %e", stability_param_);
				trcnt_ = 0;
				nstalls_ = 0;

				/** set trust region */
				setTrustRegion(stability_param_, stability_center_);
			}
			else if (nCutsAdded == 0)
			{
				nstalls_++;
			}

			message_->print(3, "\n");
		}
	}
	else
	{
		message_->print(3, "TR  dual objective %e\n", newdual);
		if (newdual >= bestdualobj_ + 1.0e-4 * (curprimobj - bestdualobj_))
			/** update dual bound */
			bestdualobj_ = newdual;

			/** update best solution */
			DSPdebugMessage("getSiPtr()->getNumCols()-nthetas_ [%d] == bestdualsol_.size() [%d]", 
				getSiPtr()->getNumCols()-nthetas_, (int) bestdualsol_.size());
			assert(getSiPtr()->getNumCols()-nthetas_==bestdualsol_.size());
			CoinCopyN(&primsol_[nthetas_], getSiPtr()->getNumCols()-nthetas_, &bestdualsol_[0]);
	}

#ifdef DSP_HAS_OOQP
	DspOsiOoqpEps * ooqp = dynamic_cast<DspOsiOoqpEps*>(osi_);
	if (ooqp)
	{
		if (ooqp->ooqp_->hasOoqpStatus_ && isSolved_)
		{
			DSPdebugMessage("bestprimobj %+e bestdualobj %+e\n", bestprimobj_, bestdualobj_);
			double epsilon = (osi_->getPrimObjValue() - newprimal + ooqp->ooqp_->getDualityGap()) / (1. + fabs(osi_->getPrimObjValue()));
			if (epsilon > 1.) epsilon = 1.;
			ooqp->ooqp_->setOoqpStatus(epsilon, -bestprimobj_, -bestdualobj_);
		}
	}
#endif

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** is solution trust region boundary? */
bool DdMasterTr::isSolutionBoundary(double eps)
{
	double maxdiff = 0.0;
	double mindiff = COIN_DBL_MAX;

	BGN_TRY_CATCH

	const double * clbd = getSiPtr()->getColLower();
	const double * cubd = getSiPtr()->getColUpper();
	int ncols = getSiPtr()->getNumCols();

	for (int j = nthetas_; j < ncols; ++j)
	{
		double diff = CoinMin(cubd[j] - primsol_[j], primsol_[j] - clbd[j]);
		if (diff > maxdiff)
			maxdiff = diff;
		if (diff < mindiff)
			mindiff = diff;
	}
	DSPdebugMessage("TR mindiff %+e\n", mindiff);

	END_TRY_CATCH(;)

	return fabs(mindiff) < eps;
}

/** add cuts */
int DdMasterTr::addCuts(
		bool possiblyDel /**< possibly delete cuts*/)
{
#define FREE_MEMORY \
	FREE_2D_ARRAY_PTR(nthetas_, aggvec); \
	FREE_ARRAY_PTR(aggrhs);

	OsiCuts cuts;
	double ** aggvec = NULL;
	double *  aggrhs = NULL;
	CoinPackedVector cutvec; /**< cut body */
	double cutrhs;           /**< cut RHS */
	int nCutsAdded; /**< number of cuts added */
	int cutidx;

	BGN_TRY_CATCH

	int ncols = getSiPtr()->getNumCols();
#ifdef DSP_DEBUG
	DSPdebugMessage("primsol_:\n");
	DspMessage::printArray(ncols, &primsol_[0]);
#endif

	/** initialize linearization error */
	linerr_ = -bestdualobj_;

	/** allocate memory for dense cut */
	aggvec = new double * [nthetas_];
	aggrhs = new double [nthetas_];
	for (int i = 0; i < nthetas_; ++i)
	{
		aggvec[i] = new double [ncols];
		CoinZeroN(aggvec[i], ncols);
		aggrhs[i] = 0.0;
	}

	/** add row cuts by looping over each scenario */
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
	{
		/** cut index */
		cutidx = s % nthetas_;

		/** calculate error and construct cut */
		linerr_ += subprimobj_[s];
		aggrhs[cutidx] += subprimobj_[s];
		for (int i = 0; i < nlambdas_; i++)
		{
			/** evaluate solution on coupling constraints (if they are Hx = d, this is (Hx - d)_i) */
			double hx_d = model_->evalLhsCouplingRowSubprob(i, s, subsolution_[s]) - model_->getRhsCouplingRow(i);
			aggvec[cutidx][nthetas_ + i] -= hx_d; /** coefficients for lambda */
			linerr_ += hx_d * stability_center_[i];
			// if (isSolved_) {
				linerr_ -= hx_d * primsol_[nthetas_ + i];
				aggrhs[cutidx] -= hx_d * primsol_[nthetas_ + i];
			// }
		}
	}

	if (model_->isDro()) {

		TssModel* tss = NULL;
		try {
			tss = dynamic_cast<TssModel*>(model_);
		} catch (const std::bad_cast &e) {
			printf("Error: Model claims to be stochastic when it is not\n");
            return DSP_RTN_ERR;
		}
		double recourse_obj;

		/** add DRO elements */
		for (int s = 0; s < model_->getNumSubproblems(); ++s) {
#ifdef DSP_DEBUG
			printf("subsolution[%d]:\n", s);
			DspMessage::printArray(tss->getNumCols(0)+tss->getNumCols(1)+1, subsolution_[s]);
#endif

			/** cut index */
			cutidx = s % nthetas_;

			// P
			recourse_obj = 0.0;
			for (int j = 0; j < tss->getNumCols(1); ++j) {
				recourse_obj += obj_reco_[s][j] * subsolution_[s][tss->getNumCols(0)+j];
				// DSPdebugMessage("obj_reco_[%d][%d] = %e, subsolution_[%d][%d] = %e\n",
				// 	s, j, obj_reco_[s][j], s, tss->getNumCols(0)+j, subsolution_[s][tss->getNumCols(0)+j]);
			}
			// recourse_obj = scen_obj->dotProduct(subsolution_[s] + tss->getNumCols(0));

			linerr_ += recourse_obj * (stability_center_[nlambdas_+nus_+s] - primsol_[ncols-nPs_+s]);
			aggvec[cutidx][ncols-nPs_+s] = -recourse_obj;
			aggrhs[cutidx] -= recourse_obj * primsol_[ncols-nPs_+s];
			DSPdebugMessage("recourse_obj[%d] = %e, primsol_[%d] = %e, aggvec[%d][%d] = %e\n",
							s, recourse_obj, ncols - nPs_ + s, primsol_[ncols - nPs_ + s], cutidx, ncols - nPs_ + s, aggvec[cutidx][ncols - nPs_ + s]);
		}
	}

	for (int s = 0; s < nthetas_; ++s)
	{
		/** construct cut */
		cutvec.clear();

		/** set it as sparse */
		aggvec[s][s] = 1.0;
		for (int j = 0; j < ncols; ++j)
		{
			if (fabs(aggvec[s][j]) > 1E-10)
				cutvec.insert(j, aggvec[s][j]);
		}

		/** cut rhs */
		cutrhs = aggrhs[s];
		if (fabs(cutrhs) < 1e-10)
			cutrhs = 0.0;
		if (model_->isDro() && fabs(cutrhs) > 1.e-10)
		{
#ifdef DSP_DEBUG
			DSPdebugMessage("cutrhs[%d] = %e\n", s, cutrhs);
#endif
			printf("Master problem may experience numerical difficulty in cut generation: (fabs(%e) >> 0.0)\n", cutrhs);
			cutrhs = 0.0;
		}
		// assert(model_->isDro() == false || cutrhs == 0.0);

		OsiRowCut * rc = new OsiRowCut;
		rc->setRow(cutvec);
		rc->setLb(-COIN_DBL_MAX);
		rc->setUb(cutrhs);
		rc->setEffectiveness(rc->violated(&primsol_[0]));
		DSPdebug(rc->print());
		DSPdebugMessage("cut effectiveness: %e\n", rc->effectiveness());

		if (rc->effectiveness() > 1.0e-6)
		{
			/** number of cuts before adding cut */
			int nCutsBefore = cuts_->sizeCuts();

			/** add cut if not duplicate */
			cuts_->insertIfNotDuplicate(*rc);

			if (nCutsBefore < cuts_->sizeCuts())
			{
				/** insertIfNotDuplicate does not set effectiveness */
				cuts_age_.push_back(0);
				possiblyDelete_.push_back(possiblyDel);
				masterobjsAtCutAdd_.push_back(osi_->getPrimObjValue());
				cuts.insert(rc);
			}
		}
		else
			FREE_PTR(rc);
	}

	nCutsAdded = cuts.sizeCuts();
	DSPdebug(cuts.printCuts());
	if (nCutsAdded > 0)
		/** apply cuts */
		getSiPtr()->applyCuts(cuts);
//	else
//		/** recruit back some cuts if no cut is generated */
//		recruiteCuts();

	END_TRY_CATCH(;)

	FREE_MEMORY

	return nCutsAdded;
#undef FREE_MEMORY
}

/** possibly delete cuts */
DSP_RTN_CODE DdMasterTr::possiblyDeleteCuts(
		double subobjval /**< sum of subproblem objective values */)
{
	BGN_TRY_CATCH

	DSP_RTN_CHECK_THROW(possiblyDeleteCutsOsi(subobjval));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** possibly delete cuts */
DSP_RTN_CODE DdMasterTr::possiblyDeleteCutsOsi(
		double subobjval /**< sum of subproblem objective values */)
{
	OsiCuts cuts;
	int nrows = model_->nonanticipativity() ? model_->getNumSubproblemCouplingCols(0) : 0;
	int ncuts = getSiPtr()->getNumRows() - nrows;
	if (ncuts == 0)
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	const double * pi = getSiPtr()->getRowPrice() + nrows;

	/** mark cuts that should not be deleted */
	for (int i = 0, i2 = 0; i < cuts_->sizeCuts(); ++i)
	{
		/** do not consider inactive cuts */
		if (cuts_age_[i] < 0) continue;
		/** consider only old cuts */
		if (cuts_age_[i] < 100)
		{
			possiblyDelete_[i] = false;
			continue;
		}
		/** aging cuts with (almost) zero Lagrangian multiplier values */
		if (fabs(pi[i2]) < 1.0e-10)
			possiblyDelete_[i] = false;
		/** do not delete cuts generated at minor iterations such that the following condition holds. */
		else if (i >= cuts_->sizeCuts() - ncuts_minor_ &&
				(getSiPtr()->getObjValue() - subobjval) > cutdel_param_ * (masterobjsAtCutAdd_[i] - subobjval))
			possiblyDelete_[i] = false;
		i2++;
	}

	/** get basis information */
	CoinWarmStartBasis * ws = NULL;
	if (getSiPtr()->getWarmStart())
		ws = dynamic_cast<CoinWarmStartBasis*>(getSiPtr()->getWarmStart()->clone());

	vector<char> aStat; /**< status of artificial variables */

	/** mark as deleted; and construct temporary cut pool to be added */
	for (int i = 0, i2 = nrows; i < cuts_->sizeCuts(); ++i)
	{
		/** do not consider inactive cuts */
		if (cuts_age_[i] < 0) continue;

		if (possiblyDelete_[i])
			cuts_age_[i] = -1;
		else
		{
			OsiRowCut * rc = cuts_->rowCutPtr(i);
			if (rc)
			{
				rc->setEffectiveness(1.0);
				cuts.insert(*rc);
				if (ws)
					aStat.push_back(ws->getArtifStatus(i2));
			}
		}

		i2++;
	}

	/** number of cuts to delete */
	int nCutsToDelete = ncuts - cuts.sizeCuts();

	/** exit if no cut to delete */
	if (nCutsToDelete == 0)
		return DSP_RTN_OK;

	/** remove all cuts from solver interface */
	removeAllCuts();

	/** apply cuts */
	getSiPtr()->applyCuts(cuts);

	if (ws)
	{
		/** create new basis */
		CoinWarmStartBasis * basis = new CoinWarmStartBasis(
				ws->getNumStructural(), ws->getNumArtificial(),
				ws->getStructuralStatus(), &aStat[0]);

		getSiPtr()->setWarmStart(basis);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** recruite cuts */
int DdMasterTr::recruiteCuts()
{
	int nRecruited = 0;
	OsiCuts cuts;

	CoinWarmStartBasis * ws = NULL;

	vector<char> aStat; /**< status of artificial variables */

	BGN_TRY_CATCH

	/** get basis information */
	getSiPtr()->setWarmStart(getSiPtr()->getWarmStart());
	ws = dynamic_cast<CoinWarmStartBasis*>(getSiPtr()->getWarmStart());

	int irow = model_->nonanticipativity() ? model_->getNumSubproblemCouplingCols(0) : 0;
	for (int i = 0; i < cuts_->sizeCuts(); ++i)
	{
		/** retrieve row cut */
		OsiRowCut * rc = cuts_->rowCutPtr(i);
		assert(rc);

		if (cuts_age_[i] >= 0)
		{
			/** add cut */
			cuts.insert(*rc);

			/** set status of artificial variable */
			aStat.push_back(ws->getArtifStatus(irow++));
		}
		else
		{
			/** set effectiveness */
			rc->setEffectiveness(rc->violated(&primsol_[0]));
			if (rc->effectiveness() > 1.0e-6)
			{
				nRecruited++;
				/** add cut */
				cuts.insert(*rc);

				/** set status of artificial variable */
				aStat.push_back(CoinWarmStartBasis::basic);

				/** other cut info */
				cuts_age_[i] = 0;
				possiblyDelete_[i] = true;
			}
		}
	}

	if (cuts.sizeCuts() > 0)
	{
		/** remove all cuts from solver interface */
		removeAllCuts();
		/** apply cuts */
		getSiPtr()->applyCuts(cuts);

		/** create new basis */
		CoinWarmStartBasis * basis = new CoinWarmStartBasis(
				ws->getNumStructural(), ws->getNumArtificial(),
				ws->getStructuralStatus(), &aStat[0]);
		getSiPtr()->setWarmStart(basis);
	}

	END_TRY_CATCH(;)

	return 0;
}

/** remove all cuts */
DSP_RTN_CODE DdMasterTr::removeAllCuts()
{
	BGN_TRY_CATCH

	int nrows = model_->nonanticipativity() ? model_->getNumSubproblemCouplingCols(0) : 0;
	int ncuts = getSiPtr()->getNumRows() - nrows;

	/** row indices to delete */
	int * rowIndices = new int [ncuts];
	CoinIotaN(rowIndices, ncuts, nrows);

	/** delete */
	getSiPtr()->deleteRows(ncuts, rowIndices);

	/** free memory */
	FREE_ARRAY_PTR(rowIndices);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** change trust region */
DSP_RTN_CODE DdMasterTr::setTrustRegion(double stability_param, double * stability_center)
{
	BGN_TRY_CATCH

	int ncols = getSiPtr()->getNumCols();
	for (int j = nthetas_; j < ncols; ++j)
	{
		double clbd = stability_center[j - nthetas_] - stability_param;
		double cubd = stability_center[j - nthetas_] + stability_param;
		if (model_->getSenseCouplingRow(j - nthetas_) == 'L')
			clbd = CoinMax((double) 0.0, clbd); /* lambda >= 0 */
		else if (model_->getSenseCouplingRow(j - nthetas_) == 'G')
			cubd = CoinMin((double) 0.0, cubd); /* lambda <= 0 */
		if (j - nthetas_ - nlambdas_ >= 0) {
			clbd = CoinMax(0.0, clbd);
			cubd = CoinMin(1.0, cubd);
		}
		getSiPtr()->setColBounds(j, clbd, cubd);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterTr::terminationTest()
{
	if (status_ == DSP_STAT_MW_STOP)
		return status_;

	int signal = status_;

	BGN_TRY_CATCH

#ifdef DSP_HAS_OOQP
	DspOsiOoqpEps * ooqp = dynamic_cast<DspOsiOoqpEps*>(osi_);
	/** is the solution suboptimal? */
	if (ooqp != NULL && ooqp->ooqp_->isSuboptimal())
		return status_;
#endif

	double absgap = getAbsApproxGap();
	double relgap = getRelApproxGap();
	DSPdebugMessage("absgap %+e relgap %+e\n", absgap, relgap);
	double gaptol = par_->getDblParam("DD/STOP_TOL");
	if (getSiPtr()->getNumIntegers() > 0)
		gaptol += par_->getDblParam("DD/SUB/GAPTOL");
	if (relgap <= gaptol) {
		signal = DSP_STAT_MW_STOP;
		status_ = DSP_STAT_OPTIMAL;
		message_->print(1, "Tr  STOP with gap tolerance %+e (%.2f%%).\n", absgap, relgap*100);
	} else if (nstalls_ >= 3 && getSiPtr()->getObjValue() < bestdualobj_) {
		signal = DSP_STAT_MW_STOP;
		status_ = DSP_STAT_STOPPED_NUMERICS;
		message_->print(1, "Tr  STOP with stalling (%d).\n", nstalls_);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return signal;
}
