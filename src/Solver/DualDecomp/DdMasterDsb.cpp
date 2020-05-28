/*
 * DdMasterDsb.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: kibaekkim
 */

#include "Solver/DualDecomp/DdMasterDsb.h"
#include "SolverInterface/DspOsi.h"

DdMasterDsb::DdMasterDsb(
		DspParams *  par,     /**< parameter pointer */
		DecModel *   model,   /**< model pointer */
		DspMessage * message /**< message pointer */) :
DdMaster(par, model, message),
prox_(NULL),
phi_t_(10.0),
phi_l_(1000.0),
alpha_t_(0.1),
alpha_l_(0.1),
eps_opt_(1.0e-6),
eps_e_(1.0e-6),
eps_g_(1.0e-6),
tau_min_(1.0e-6),
tau_(1.0),
upperbound_(1.0e+20),
valueAtProx_(-1.0e+20),
level_(0.0),
gg_(0.0),
isSolved_(false),
nrows_(0),
ncols_(0),
nzcnt_(0),
modelObjval_(0.0),
bb_(0.0) {}

DdMasterDsb::~DdMasterDsb()
{
	FREE_ARRAY_PTR(prox_);
	FREE_ARRAY_PTR(lambda_);
	for (unsigned int i = 0; i < lhsCouplingRows_.size(); ++i)
		FREE_ARRAY_PTR(lhsCouplingRows_[i]);
}

/** initialize */
DSP_RTN_CODE DdMasterDsb::init()
{
	BGN_TRY_CATCH

	DdMaster::init();

	/** create problem */
	createProblem();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterDsb::solve()
{
	BGN_TRY_CATCH

	bool resolve = false;
	int tmpcnt = 0;
	double cputime  = CoinCpuTime();
	double walltime = CoinGetTimeOfDay();
	int status = 0;
	double objval = 0.0;

	while (1)
	{
		resolve = false;

		/** solve */
		si_->resolve();

		/** mark as solved */
		isSolved_ = true;

		if (si_->isProvenOptimal()) {
			status = DSP_STAT_OPTIMAL;
			/** TODO realloc */
			FREE_ARRAY_PTR(primsol_);
			primsol_ = new double [si_->getNumCols()];

			/** copy solution */
			CoinCopyN(si_->getColSolution(), si_->getNumCols(), primsol_);
#ifdef DSP_DEBUG
			DSPdebugMessage("Master solution:");
			for (int j = 0; j < ncols; ++j)
			{
				if (j % 5 == 0) printf("\n\t");
				printf("%e ", primsol_[j]);
			}
			printf("\n");
#endif

			/** construct dual solution:
			 *    General:
			 *      \lambda = prox_k + tau_k * (b u + \sum_j \sum_s B_s x_s^j y_s^j)
			 *    Stochastic:
			 *      \lambda = prox_k + tau_k * (w + \sum_j x_s^j y_s^j) */
			if (model_->nonanticipativity())
			{
				for (int i = 0; i < model_->getNumCouplingRows(); ++i)
					lambda_[i] = prox_[i] + tau_ * primsol_[2 + (i % model_->getNumSubproblemCouplingCols(0))];

				for (unsigned int j = 0, pos = 0; j < lhsCouplingRows_.size(); ++j)
				{
					if (activeCols_[j] == false) continue;
					int sid = dynids_[j];
					for (int i = 0; i < model_->getNumSubproblemCouplingCols(0); ++i)
						lambda_[model_->getNumSubproblemCouplingCols(0) * sid + i] += tau_ * primsol_[si_->getNumCols() + pos] * lhsCouplingRows_[j][model_->getNumSubproblemCouplingCols(0) * sid + i];
					pos++;
				}
			}
			else
			{
				for (int i = 0; i < model_->getNumCouplingRows(); ++i)
				{
					double yBx = 0.0;
					for (unsigned int j = 0, pos = 0; j < lhsCouplingRows_.size(); ++j)
						if (activeCols_[j])
						{
							yBx += primsol_[si_->getNumCols() + pos] * lhsCouplingRows_[j][i];
							pos++;
						}
					lambda_[i] = prox_[i] + tau_ * (primsol_[0] * model_->getRhsCouplingRow(i) + yBx);
				}
			}
#ifdef DSP_DEBUG
			DSPdebugMessage("Lambda:");
			for (int i = 0; i < ncoupling_; ++i)
			{
				if (i % 5 == 0) printf("\n\t");
				printf("%e ", lambda_[i]);
			}
			printf("\n");
#endif

			objval = si_->getObjValue();
			for (int i = 0; i < model_->getNumCouplingRows(); ++i)
				objval += pow(lambda_[i] - prox_[i], 2.0) / (2.0 * tau_);
			DSPdebugMessage("Original objective of the bundle problem: %e\n", objval);
			DSPdebugMessage("Manually calculated objective value without proximal term: %e\n", objval);
			primobj_ = objval;

			/** update phi_t */
			phi_t_ = primobj_ - valueAtProx_;

			/** 2-norm of subgradient */
			gg_ = (primobj_ - objval) * 2 * tau_;
		} else {
			/** get status */
			if (si_->isIterationLimitReached())
				status = DSP_STAT_LIM_ITERorTIME;
			else if (si_->isProvenPrimalInfeasible())
				status = DSP_STAT_PRIM_INFEASIBLE;
			else if (si_->isProvenDualInfeasible())
				status = DSP_STAT_DUAL_INFEASIBLE;
			/** zero solution */
			CoinZeroN(primsol_, si_->getNumCols());

			upperbound_ = min(level_, upperbound_);
			phi_l_ = max(eps_opt_, (1 - alpha_l_) * (upperbound_ - valueAtProx_));
			level_ = valueAtProx_ + phi_l_;
			message_->print(3, "Updates upperbound %e phi_l %e level %e\n", upperbound_, phi_l_, level_);
			message_->print(2, "  Warning: solution status %d\n", status);
			message_->print(2, "           UB %e phi_l %e level %e\n", upperbound_, phi_l_, level_);

			/** update linear objective coefficient */
			applyLevelChange();

			/** mark as resolve */
			if (upperbound_ - valueAtProx_ >= eps_opt_)
				resolve = true;
		}

		if (!resolve || tmpcnt > 5)
			break;
	}

	/** update statistics */
	s_statuses_.push_back(status);
	s_primobjs_.push_back(objval);
	s_dualobjs_.push_back(objval);
	double * s_primsol = new double [si_->getNumCols()];
	CoinCopyN(primsol_, si_->getNumCols(), s_primsol);
	s_primsols_.push_back(s_primsol);
	s_primsol = NULL;
	s_cputimes_.push_back(CoinCpuTime() - cputime);
	s_walltimes_.push_back(CoinGetTimeOfDay() - walltime);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterDsb::createProblem()
{
#define FREE_MEMORY            \
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
		FREE_ARRAY_PTR(irowQ); \
		FREE_ARRAY_PTR(jcolQ); \
		FREE_ARRAY_PTR(dQ);    \
		FREE_PTR(mat);

	char * ctype  = NULL; /**< column type */
	double * obj  = NULL; /**< linear objective function coefficient */
	double * clbd = NULL; /**< column lower bound */
	double * cubd = NULL; /**< column upper bound */
	double * rlbd = NULL; /**< row lower bound */
	double * rubd = NULL; /**< row upper bound */

	/** linear constraint matrix */
	CoinBigIndex * bgn = NULL; /**< matrix start index */
	int * len = NULL;
	int * ind = NULL;
	double * elem = NULL;
	CoinPackedMatrix * mat = NULL;

	/** Hessian matrix */
	int nnzQ = 0;       /**< number of nonzeros */
	int * irowQ = NULL; /**< row indicies of Hessian */
	int * jcolQ = NULL; /**< column indicies of Hessian */
	double * dQ = NULL; /**< elements of Hessian */

	BGN_TRY_CATCH

	/** LP dimension */
	nrows_ = model_->getNumSubproblems() + 1;
	ncols_ = 2;
	nzcnt_  = model_->getNumSubproblems() + 2;
	if (model_->nonanticipativity())
		ncols_ += model_->getNumSubproblemCouplingCols(0);
	else
	{
		for (int i = 0; i < model_->getNumCouplingRows(); ++i)
			if (model_->getSenseCouplingRow(i) != 'E')
			{
				nrows_++;
				nzcnt_++;
			}
	}

	/** allocate memory */
	ctype = new char [ncols_];
	clbd = new double[ncols_];
	cubd = new double[ncols_];
	obj  = new double[ncols_];
	rlbd = new double[nrows_];
	rubd = new double[nrows_];
	bgn  = new CoinBigIndex[nrows_ + 1];
	len  = new int[nrows_];
	ind  = new int[nzcnt_];
	elem = new double[nzcnt_];

	if (model_->nonanticipativity())
	{
		irowQ = new int [model_->getNumSubproblemCouplingCols(0)];
		jcolQ = new int [model_->getNumSubproblemCouplingCols(0)];
		dQ    = new double [model_->getNumSubproblemCouplingCols(0)];
	}
	else
	{
		irowQ = new int [1];
		jcolQ = new int [1];
		dQ    = new double [1];
	}

	/** all continuous variables */
	CoinFillN(ctype, ncols_, 'C');

	/** proximal point */
	prox_ = new double [model_->getNumCouplingRows()];
	CoinZeroN(prox_, model_->getNumCouplingRows());

	/** dual variable */
	lambda_  = new double [model_->getNumCouplingRows()];
	CoinZeroN(lambda_, model_->getNumCouplingRows());

	/** linear objective coefficient */
	obj[0] = 0.0;
	obj[1] = -level_;
	if (model_->nonanticipativity())
	{
		CoinZeroN(obj + 2, model_->getNumSubproblemCouplingCols(0));
		for (int s = 0; s < model_->getNumSubproblems(); ++s)
		{
			for (int j = 0; j < model_->getNumSubproblemCouplingCols(0); ++j)
				obj[2 + j] += prox_[model_->getNumSubproblemCouplingCols(0) * s + j];
		}
	}

	/** bounds */
	if (model_->nonanticipativity())
	{
		CoinFillN(clbd, 2, 0.0);
		CoinFillN(clbd + 2, ncols_ - 2, -COIN_DBL_MAX);
	}
	else
		CoinFillN(clbd, ncols_, 0.0);
	CoinFillN(cubd, ncols_, COIN_DBL_MAX);

	/** row bounds */
	CoinFillN(rlbd, model_->getNumSubproblems(), 0.0);
	CoinFillN(rubd, model_->getNumSubproblems(), 0.0);
	CoinFillN(rlbd + model_->getNumSubproblems(), 1, 1.0);
	CoinFillN(rubd + model_->getNumSubproblems(), 1, 1.0);
	if (model_->nonanticipativity() == false)
	{
		for (int i = 0, ii = 0; i < model_->getNumCouplingRows(); ++i)
		{
			switch(model_->getSenseCouplingRow(i))
			{
			case 'L':
				rlbd[model_->getNumSubproblems()+1+ii] = -COIN_DBL_MAX;
				rubd[model_->getNumSubproblems()+1+ii] = -prox_[i] / tau_;
				ii++;
				break;
			case 'G':
				rlbd[model_->getNumSubproblems()+1+ii] = -prox_[i] / tau_;
				rubd[model_->getNumSubproblems()+1+ii] = COIN_DBL_MAX;
				ii++;
				break;
			default:
				break;
			}
		}
	}

	/** for constraints */
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
	{
		bgn[i] = i;
		ind[i] = 0;
		elem[i] = 1.0;
		len[i] = 1;
	}
	bgn[model_->getNumSubproblems()] = model_->getNumSubproblems();
	ind[model_->getNumSubproblems()] = 0;
	elem[model_->getNumSubproblems()] = 1.0;
	ind[model_->getNumSubproblems()+1] = 1;
	elem[model_->getNumSubproblems()+1] = -1.0;
	len[model_->getNumSubproblems()] = 2;
	if (model_->nonanticipativity() == false)
	{
		for (int i = 0, ii = 0; i < model_->getNumCouplingRows(); ++i)
		{
			if (model_->getSenseCouplingRow(i) != 'E')
			{
				bgn[model_->getNumSubproblems()+1+ii] = model_->getNumSubproblems() + 2 + ii;
				ind[model_->getNumSubproblems()+2+ii] = 0;
				elem[model_->getNumSubproblems()+2+ii] = model_->getRhsCouplingRow(i);
				len[model_->getNumSubproblems()+1+ii] = 1;
				ii++;
			}
		}
	}

	/** constraint matrix */
	mat = new CoinPackedMatrix(false, ncols_, nrows_, nzcnt_, elem, ind, bgn, len);
	DSPdebug(mat->verifyMtx(4));

	/** create solver interface */
	switch (par_->getIntParam("DD/MASTER_ALGO"))
	{
	case DSBM:
		si_ = new DspOsiOoqp();
		si_->setLogLevel(par_->getIntParam("LOG_LEVEL"));
		break;
	default:
		break;
	}

	/** copy problem data */
	si_->loadProblem(*mat, clbd, cubd, obj, ctype, rlbd, rubd);

	if (model_->nonanticipativity())
	{
		for (int j = 0; j < model_->getNumSubproblemCouplingCols(0); ++j)
		{
			irowQ[j] = 2 + j;
			jcolQ[j] = 2 + j;
			dQ[j] = tau_ * model_->getNumSubproblems();
		}
		nnzQ = model_->getNumSubproblemCouplingCols(0);
	}
	else
	{
		/** calculate b^T b */
		bb_ = 0.0;
		for (int i = 0; i < model_->getNumCouplingRows(); ++i)
			bb_ += model_->getRhsCouplingRow(i) * model_->getRhsCouplingRow(i);

		/** set hessian */
		irowQ[0] = 0; jcolQ[0] = 0; dQ[0] = tau_ * model_->getNumSubproblems() * bb_;
		if (fabs(dQ[0]) > 1.0e-10)
			nnzQ = 1;
	}
	SolverInterfaceOoqp * ooqp = dynamic_cast<SolverInterfaceOoqp*>(si_);
	if (ooqp)
		ooqp->setHessian(nnzQ, irowQ, jcolQ, dQ);

	/** allocate memory for solution */
	primsol_ = new double [ncols_];
	CoinZeroN(primsol_, ncols_);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY;

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMasterDsb::updateProblem()
{
	BGN_TRY_CATCH

	/** update upper bound */
	if(bestprimobj_ < upperbound_)
		upperbound_ = bestprimobj_;

	/** current objective value */
	double curObjval = 0.0;
	for (int i = 0; i < model_->getNumSubproblems(); ++i)
		curObjval += subdualobj_[i];
	message_->print(3, "Current objective value: %f\n", curObjval);

	/** update LHS of coupling rows
	 * and dynamic objective coefficients */
	for (int s = 0; s < model_->getNumSubproblems(); ++s)
	{
		double obj = subprimobj_[s];
		double * lhs = new double [model_->getNumCouplingRows()];
		for (int i = 0; i < model_->getNumCouplingRows(); ++i)
		{
			lhs[i] = model_->evalLhsCouplingRowSubprob(i, s, subsolution_[s]);
			obj += (prox_[i] - lambda_[i]) * lhs[i];
		}
		dynObjs_.push_back(obj);
		lhsCouplingRows_.push_back(lhs);
		dynids_.push_back(s);
	}

	/** ascent step acceptance test */
	if (curObjval >= valueAtProx_ + alpha_t_ * phi_t_)
	{
		valueAtProx_ = curObjval;
		bestdualobj_ = valueAtProx_;

		/** update objective coefficient of dynamic columns */
		updateDynObjs();

		if (isSolved_)
		{
			/** update proximal point and its value */
			CoinCopyN(lambda_, model_->getNumCouplingRows(), prox_);

			/** update tau_ */
			tau_ *= primsol_[0];

			/** TODO update RHS */
		}
		/** update phi_l_ */
		phi_l_ = min(phi_l_, (1-alpha_l_) * (upperbound_ - valueAtProx_));
		message_->print(3, "Ascent step: tau %f phi_l %f\n", tau_, phi_l_);
	}
	else
	{
		/** TODO: may change tau_ */
		tau_ = max(1.0e-4, tau_ * phi_l_ / phi_t_);
		DSPdebugMessage("Null step: updates tau_ %e\n", tau_);

		if (isSolved_)
		{
			/** update phi_l_ */
			if (primsol_[0] > 1)
			{
				phi_l_ = max(phi_l_, phi_l_ * alpha_l_);
				message_->print(3, "Null step: updates phi_l %f\n", phi_l_);
			}
		}
	}

	/** manage bundle */
	manageBundle(subprimobj_);

	/** update level_ */
	level_ = valueAtProx_ + phi_l_;
	message_->print(1, "  Level updated %f\n", level_);

	/** update objective function */
	applyLevelChange();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterDsb::manageBundle(double* objvals)
{
	BGN_TRY_CATCH

	int posi, posj;
	vector<int> irowQ;
	vector<int> jcolQ;
	vector<double> dQ;
	SolverInterfaceOoqp * ooqp = dynamic_cast<SolverInterfaceOoqp*>(si_);

	/** TODO: delete some columns */
	ooqp->clearCols();
	activeCols_.clear();

	/** number of columns added */
	int ncolsAdded = activeCols_.size();

	/** add new columns */
	for (int j = ncolsAdded; j < dynObjs_.size(); ++j)
	{
		int sid = dynids_[j];
		vector<int> indices;
		vector<double> vals;

		indices.push_back(sid);
		vals.push_back(-1.0);

		if (model_->nonanticipativity() == false)
		{
			for (int i = 0; i < model_->getNumCouplingRows(); ++i)
			{
				if (model_->getSenseCouplingRow(i) != 'E')
				{
					if (fabs(lhsCouplingRows_[j][i]) < 1.0e-10)
						continue;
					indices.push_back(model_->getNumSubproblems() + 1 + i);
					vals.push_back(lhsCouplingRows_[j][i]);
				}
			}
		}

		ooqp->addCol(indices.size(), &indices[0], &vals[0], 0.0, COIN_DBL_MAX, dynObjs_[j]);

		/** set as active */
		activeCols_.push_back(true);
	}

	/** re-create Hessian */
	if (model_->nonanticipativity())
	{
		for (int j = 0; j < model_->getNumSubproblemCouplingCols(0); ++j)
		{
			irowQ.push_back(2 + j);
			jcolQ.push_back(2 + j);
			dQ.push_back(tau_ * model_->getNumSubproblems());
		}
	}
	else if (fabs(bb_) > 1.0e-10)
	{
		irowQ.push_back(0);
		jcolQ.push_back(0);
		dQ.push_back(tau_ * model_->getNumSubproblems() * bb_);
	}

	posi = si_->getNumCols();

	for (unsigned int j = 0; j < activeCols_.size(); ++j)
	{
		if (activeCols_[j] == false) continue;

		int sid = dynids_[j];

		if (model_->nonanticipativity())
		{
			/** This sets x_s^j in the Hessian matrix. */
			for (int i = 0; i < model_->getNumSubproblemCouplingCols(0); ++i)
			{
				if (fabs(lhsCouplingRows_[j][model_->getNumSubproblemCouplingCols(0) * sid + i]) < 1.0e-10)
					continue;
				irowQ.push_back(posi);
				jcolQ.push_back(2 + i);
				dQ.push_back(tau_ * lhsCouplingRows_[j][model_->getNumSubproblemCouplingCols(0) * sid + i]);
			}
		}
		else
		{
			/** This sets (b^T B_s x_s^j) in the Hessian matrix. */
			double val = 0.0;
			for (int i = 0; i < model_->getNumCouplingRows(); ++i)
				val += model_->getRhsCouplingRow(i) * lhsCouplingRows_[j][i];

			if (fabs(val) > 1.0e-10)
			{
				irowQ.push_back(posi);
				jcolQ.push_back(0);
				dQ.push_back(tau_ * val);
			}
		}
		posj = si_->getNumCols();

		for (unsigned int k = 0; k < activeCols_.size(); ++k)
		{
			if (k > j) break;
			if (activeCols_[k] == false) continue;

			if (dynids_[k] != sid)
			{
				posj++;
				continue;
			}

			double val = 0.0;
			if (model_->nonanticipativity())
			{
				int sid2 = dynids_[k];
				for (int i = 0; i < model_->getNumSubproblemCouplingCols(0); ++i)
					val += lhsCouplingRows_[j][model_->getNumSubproblemCouplingCols(0) * sid + i] * lhsCouplingRows_[k][model_->getNumSubproblemCouplingCols(0) * sid2 + i];
			}
			else
			{
				for (int i = 0; i < model_->getNumCouplingRows(); ++i)
					val += lhsCouplingRows_[j][i] * lhsCouplingRows_[k][i];
			}
			//DSPdebugMessage("irowQ %d jcolQ %d dQ %f\n", posi, posj, tau_*val);

			if (fabs(val) > 1.0e-10)
			{
				irowQ.push_back(posi);
				jcolQ.push_back(posj);
				dQ.push_back(tau_ * val);
			}

			posj++;
		}

		posi++;
	}

#ifdef DSP_DEBUG
		DSPdebugMessage("Hessian Q:");
		for (unsigned i = 0; i < irowQ.size(); ++i)
		{
			if (i % 5 == 0) printf("\n\t");
			printf("(%d,%d)[%e] ", irowQ[i], jcolQ[i], dQ[i]);
		}
		printf("\n");
#endif
	/** update Hessian */
	ooqp->setHessian(irowQ.size(), &irowQ[0], &jcolQ[0], &dQ[0]);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DdMasterDsb::applyLevelChange()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(obj);

	double * obj = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	obj = new double [si_->getNumCols()];

	CoinZeroN(obj, si_->getNumCols());
	obj[1] = -level_;

	si_->setObjective(obj);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY;

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DdMasterDsb::updateDynObjs()
{
	BGN_TRY_CATCH

	int nobjs = dynObjs_.size();

	for (int k = 0; k < nobjs; ++k)
	{
		for (int i = 0; i < model_->getNumCouplingRows(); ++i)
			dynObjs_[k] += lhsCouplingRows_[k][i] * (lambda_[i] - prox_[i]);
		DSPdebugMessage("dynObjs_[%d] %e\n", k, dynObjs_[k]);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
