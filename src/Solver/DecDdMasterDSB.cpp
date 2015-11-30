/*
 * DecDdMasterDSB.cpp
 *
 *  Created on: Aug 12, 2015
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "Utility/StoMessage.h"
#include "Utility/DspParams.h"
#include <Solver/DecDdMasterDSB.h>

DecDdMasterDSB::~DecDdMasterDSB()
{
	FREE_PTR(si_);
	model_ = NULL;
	FREE_ARRAY_PTR(prox_);
	FREE_ARRAY_PTR(lambda_);
	for (unsigned int i = 0; i < lhsCouplingRows_.size(); ++i)
		FREE_ARRAY_PTR(lhsCouplingRows_[i]);
}

STO_RTN_CODE DecDdMasterDSB::createProblem(DecModel * model)
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

	/** retrieve model information */
	model_     = model;
	nsubprobs_ = model_->getNumSubproblems();
	ncoupling_ = model_->getNumCouplingRows();
	ncols0_    = model_->getNumSubproblemCouplingCols(0);

	/** LP dimension */
	nrows_ = nsubprobs_ + 1;
	ncols_ = 2;
	nzcnt_  = nsubprobs_ + 2;
	if (model_->nonanticipativity())
		ncols_ += ncols0_;
	else
	{
		for (int i = 0; i < ncoupling_; ++i)
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
		irowQ = new int [ncols0_];
		jcolQ = new int [ncols0_];
		dQ    = new double [ncols0_];
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
	prox_ = new double [ncoupling_];
	CoinZeroN(prox_, ncoupling_);

	/** dual variable */
	lambda_  = new double [ncoupling_];
	CoinZeroN(lambda_, ncoupling_);

	/** linear objective coefficient */
	obj[0] = 0.0;
	obj[1] = -level_;
	if (model_->nonanticipativity())
	{
		CoinZeroN(obj + 2, ncols0_);
		for (int s = 0; s < nsubprobs_; ++s)
		{
			for (int j = 0; j < ncols0_; ++j)
				obj[2 + j] += prox_[ncols0_ * s + j];
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
	CoinFillN(rlbd, nsubprobs_, 0.0);
	CoinFillN(rubd, nsubprobs_, 0.0);
	CoinFillN(rlbd + nsubprobs_, 1, 1.0);
	CoinFillN(rubd + nsubprobs_, 1, 1.0);
	if (model_->nonanticipativity() == false)
	{
		for (int i = 0, ii = 0; i < ncoupling_; ++i)
		{
			switch(model_->getSenseCouplingRow(i))
			{
			case 'L':
				rlbd[nsubprobs_+1+ii] = -COIN_DBL_MAX;
				rubd[nsubprobs_+1+ii] = -prox_[i] / tau_;
				ii++;
				break;
			case 'G':
				rlbd[nsubprobs_+1+ii] = -prox_[i] / tau_;
				rubd[nsubprobs_+1+ii] = COIN_DBL_MAX;
				ii++;
				break;
			default:
				break;
			}
		}
	}

	/** for constraints */
	for (int i = 0; i < nsubprobs_; ++i)
	{
		bgn[i] = i;
		ind[i] = 0;
		elem[i] = 1.0;
		len[i] = 1;
	}
	bgn[nsubprobs_] = nsubprobs_;
	ind[nsubprobs_] = 0;
	elem[nsubprobs_] = 1.0;
	ind[nsubprobs_+1] = 1;
	elem[nsubprobs_+1] = -1.0;
	len[nsubprobs_] = 2;
	if (model_->nonanticipativity() == false)
	{
		for (int i = 0, ii = 0; i < ncoupling_; ++i)
		{
			if (model_->getSenseCouplingRow(i) != 'E')
			{
				bgn[nsubprobs_+1+ii] = nsubprobs_ + 2 + ii;
				ind[nsubprobs_+2+ii] = 0;
				elem[nsubprobs_+2+ii] = model_->getRhsCouplingRow(i);
				len[nsubprobs_+1+ii] = 1;
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
		si_ = new SolverInterfaceOoqp(par_);
		si_->setPrintLevel(par_->getIntParam("LOG_LEVEL"));
		break;
	default:
		break;
	}

	/** copy problem data */
	si_->loadProblem(mat, clbd, cubd, obj, ctype, rlbd, rubd);

	if (model_->nonanticipativity())
	{
		for (int j = 0; j < ncols0_; ++j)
		{
			irowQ[j] = 2 + j;
			jcolQ[j] = 2 + j;
			dQ[j] = tau_ * nsubprobs_;
		}
		nnzQ = ncols0_;
	}
	else
	{
		/** calculate b^T b */
		bb_ = 0.0;
		for (int i = 0; i < ncoupling_; ++i)
			bb_ += model_->getRhsCouplingRow(i) * model_->getRhsCouplingRow(i);

		/** set hessian */
		irowQ[0] = 0; jcolQ[0] = 0; dQ[0] = tau_ * nsubprobs_ * bb_;
		if (fabs(dQ[0]) > 1.0e-10)
			nnzQ = 1;
	}
	si_->setHessian(nnzQ, irowQ, jcolQ, dQ);

	/** allocate memory for solution */
	solution_ = new double [ncols_];
	CoinZeroN(solution_, ncols_);

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** update problem: may update dual bound */
STO_RTN_CODE DecDdMasterDSB::updateProblem(
		double primal_bound,     /**< primal bound of the original problem */
		double & dual_bound,     /**< dual bound of the original problem */
		double * primal_objvals, /**< objective values of subproblems */
		double * dual_objvals,   /**< objective values of subproblems */
		double ** solution       /**< subproblem solutions */)
{
	BGN_TRY_CATCH

	/** update upper bound */
	if(primal_bound < upperbound_)
		upperbound_ = primal_bound;

	/** current objective value */
	double curObjval = 0.0;
	for (int i = 0; i < nsubprobs_; ++i)
		curObjval += dual_objvals[i];
	DSPdebugMessage("Current objective value: %f\n", curObjval);

	/** update LHS of coupling rows
	 * and dynamic objective coefficients */
	for (int s = 0; s < nsubprobs_; ++s)
	{
		double obj = primal_objvals[s];
		double * lhs = new double [ncoupling_];
		for (int i = 0; i < ncoupling_; ++i)
		{
			lhs[i] = model_->evalLhsCouplingRowSubprob(i, s, solution[s]);
			obj += (prox_[i] - lambda_[i]) * lhs[i];
		}
		//DSPdebugMessage("s %d obj %e\n", s, obj);
		dynObjs_.push_back(obj);
		lhsCouplingRows_.push_back(lhs);
		dynids_.push_back(s);
	}

	/** ascent step acceptance test */
	if (curObjval >= valueAtProx_ + alpha_t_ * phi_t_)
	{
		valueAtProx_ = curObjval;
		dual_bound = valueAtProx_;

		/** update objective coefficient of dynamic columns */
		updateDynObjs();

		if (isSolved_)
		{
			/** update proximal point and its value */
			CoinCopyN(lambda_, ncoupling_, prox_);

			/** update tau_ */
			tau_ *= solution_[0];

			/** TODO update RHS */
		}
		/** update phi_l_ */
		phi_l_ = min(phi_l_, (1-alpha_l_) * (upperbound_ - valueAtProx_));
		DSPdebugMessage("Ascent step: tau %f phi_l %f\n", tau_, phi_l_);
	}
	else
	{
		/** TODO: may change tau_ */
		tau_ = max(1.0e-4, tau_ * phi_l_ / phi_t_);
		DSPdebugMessage("Null step: updates tau_ %e\n", tau_);

		if (isSolved_)
		{
			/** update phi_l_ */
			if (solution_[0] > 1)
			{
				phi_l_ = max(phi_l_, phi_l_ * alpha_l_);
				DSPdebugMessage("Null step: updates phi_l %f\n", phi_l_);
			}
		}
	}

	/** manage bundle */
	manageBundle(primal_objvals);

	/** update level_ */
	level_ = valueAtProx_ + phi_l_;
	if (par_->getIntParam("LOG_LEVEL") > 1)
		printf("  Level updated %f\n", level_);

	/** update objective function */
	applyLevelChange();

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** apply level change */
STO_RTN_CODE DecDdMasterDSB::applyLevelChange()
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(obj);

	double * obj = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	obj = new double [si_->getNumCols()];

	CoinZeroN(obj, si_->getNumCols());
	obj[1] = -level_;

	si_->setObjCoef(obj);

	END_TRY_CATCH_RTN(FREE_MEMORY,STO_RTN_ERR)

	FREE_MEMORY;

	return STO_RTN_OK;
#undef FREE_MEMORY
}

/** update dynamic objective coefficients */
STO_RTN_CODE DecDdMasterDSB::updateDynObjs()
{
	BGN_TRY_CATCH

	int nobjs = dynObjs_.size();

	for (int k = 0; k < nobjs; ++k)
	{
		for (int i = 0; i < ncoupling_; ++i)
			dynObjs_[k] += lhsCouplingRows_[k][i] * (lambda_[i] - prox_[i]);
		DSPdebugMessage("dynObjs_[%d] %e\n", k, dynObjs_[k]);
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** solve problem */
STO_RTN_CODE DecDdMasterDSB::solve()
{
	BGN_TRY_CATCH

	bool resolve = false;
	int tmpcnt = 0;

	while (1)
	{
		resolve = false;

		/** solve */
		si_->solve();

		/** mark as solved */
		isSolved_ = true;

		switch(si_->getStatus())
		{
		case STO_STAT_OPTIMAL:
		{
			int ncols = si_->getNumCols() + si_->dynCols_->size();

			/** TODO realloc */
			FREE_ARRAY_PTR(solution_);
			solution_ = new double [ncols];

			/** copy solution */
			CoinCopyN(si_->getSolution(), ncols, solution_);
#ifdef DSP_DEBUG
			DSPdebugMessage("Master solution:");
			for (int j = 0; j < ncols; ++j)
			{
				if (j % 5 == 0) printf("\n\t");
				printf("%e ", solution_[j]);
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
				for (int i = 0; i < ncoupling_; ++i)
					lambda_[i] = prox_[i] + tau_ * solution_[2 + (i % ncols0_)];

				for (unsigned int j = 0, pos = 0; j < lhsCouplingRows_.size(); ++j)
				{
					if (activeCols_[j] == false) continue;
					int sid = dynids_[j];
					for (int i = 0; i < ncols0_; ++i)
						lambda_[ncols0_ * sid + i] += tau_ * solution_[si_->getNumCols() + pos] * lhsCouplingRows_[j][ncols0_ * sid + i];
					pos++;
				}
			}
			else
			{
				for (int i = 0; i < ncoupling_; ++i)
				{
					double yBx = 0.0;
					for (unsigned int j = 0, pos = 0; j < lhsCouplingRows_.size(); ++j)
						if (activeCols_[j])
						{
							yBx += solution_[si_->getNumCols() + pos] * lhsCouplingRows_[j][i];
							pos++;
						}
					lambda_[i] = prox_[i] + tau_ * (solution_[0] * model_->getRhsCouplingRow(i) + yBx);
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

			double objval = si_->getPrimalBound();
			for (int i = 0; i < ncoupling_; ++i)
				objval += pow(lambda_[i] - prox_[i], 2.0) / (2.0 * tau_);
			DSPdebugMessage("Original objective of the bundle problem: %e\n", si_->getPrimalBound());
			DSPdebugMessage("Manually calculated objective value without proximal term: %e\n", objval);
			objval_ = objval;

			/** objective value without proximal term */
//			objval_ = 0.0;
//			for (int j = 0; j < si_->getNumCols(); ++j)
//				objval_ += si_->getObjCoef()[j] * solution_[j];
//			for (int j = 0; j < si_->dynCols_->size(); ++j)
//				objval_ += si_->dynCols_->objs_[j] * solution_[si_->getNumCols() + j];
			objvals_.push_back(objval_);

			for (int j = 0; j < nsubprobs_ + 1; ++j)
				DSPdebugMessage("j %d y %e\n", j, si_->y()[j]);
			/** update phi_t */
			phi_t_ = objval_ - valueAtProx_;
			DSPdebugMessage("phi_t %e\n", phi_t_);

			/** 2-norm of subgradient */
			gg_ = (objval_ - si_->getPrimalBound()) * 2 * tau_;

			break;
		}
		default:
		{
			upperbound_ = min(level_, upperbound_);
			phi_l_ = max(eps_opt_, (1 - alpha_l_) * (upperbound_ - valueAtProx_));
			level_ = valueAtProx_ + phi_l_;
			//level_ = -108390;
			DSPdebugMessage("Updates upperbound %e phi_l %e level %e\n", upperbound_, phi_l_, level_);
			if (par_->getIntParam("LOG_LEVEL") > 2)
			{
				printf("  Warning: solution status %d\n", si_->getStatus());
				printf("           UB %e phi_l %e level %e\n", upperbound_, phi_l_, level_);
			}

			/** update linear objective coefficient */
			applyLevelChange();

			/** mark as resolve */
			if (upperbound_ - valueAtProx_ >= eps_opt_)
				resolve = true;
//			tmpcnt++;

			break;
		}
		}

		if (!resolve || tmpcnt > 5)
			break;
	}

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** manage bundle */
STO_RTN_CODE DecDdMasterDSB::manageBundle(
		double * objvals /**< subproblem objective values */)
{
	BGN_TRY_CATCH

	int posi, posj;
	vector<int> irowQ;
	vector<int> jcolQ;
	vector<double> dQ;

	/** TODO: delete some columns */
	si_->clearCols();
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
			for (int i = 0; i < ncoupling_; ++i)
			{
				if (model_->getSenseCouplingRow(i) != 'E')
				{
					if (fabs(lhsCouplingRows_[j][i]) < 1.0e-10)
						continue;
					indices.push_back(nsubprobs_ + 1 + i);
					vals.push_back(lhsCouplingRows_[j][i]);
				}
			}
		}
#ifdef DSP_DEBUG
		DSPdebugMessage("Column %d (obj %e):", j, dynObjs_[j]);
		for (unsigned i = 0; i < indices.size(); ++i)
		{
			if (i % 5 == 0) printf("\n\t");
			printf("%d [%e] ", indices[i], vals[i]);
		}
		printf("\n");
#endif

		si_->addCol(indices.size(), &indices[0], &vals[0], 0.0, COIN_DBL_MAX, dynObjs_[j]);

		/** set as active */
		activeCols_.push_back(true);
	}

	/** re-create Hessian */
	if (model_->nonanticipativity())
	{
		for (int j = 0; j < ncols0_; ++j)
		{
			irowQ.push_back(2 + j);
			jcolQ.push_back(2 + j);
			dQ.push_back(tau_ * nsubprobs_);
		}
	}
	else if (fabs(bb_) > 1.0e-10)
	{
		irowQ.push_back(0);
		jcolQ.push_back(0);
		dQ.push_back(tau_ * nsubprobs_ * bb_);
	}

	posi = si_->getNumCols();

	for (unsigned int j = 0; j < activeCols_.size(); ++j)
	{
		if (activeCols_[j] == false) continue;

		int sid = dynids_[j];

		if (model_->nonanticipativity())
		{
			/** This sets x_s^j in the Hessian matrix. */
			for (int i = 0; i < ncols0_; ++i)
			{
				if (fabs(lhsCouplingRows_[j][ncols0_ * sid + i]) < 1.0e-10)
					continue;
				irowQ.push_back(posi);
				jcolQ.push_back(2 + i);
				dQ.push_back(tau_ * lhsCouplingRows_[j][ncols0_ * sid + i]);
			}
		}
		else
		{
			/** This sets (b^T B_s x_s^j) in the Hessian matrix. */
			double val = 0.0;
			for (int i = 0; i < ncoupling_; ++i)
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
				for (int i = 0; i < ncols0_; ++i)
					val += lhsCouplingRows_[j][ncols0_ * sid + i] * lhsCouplingRows_[k][ncols0_ * sid2 + i];
			}
			else
			{
				for (int i = 0; i < ncoupling_; ++i)
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
	si_->setHessian(irowQ.size(), &irowQ[0], &jcolQ[0], &dQ[0]);

	END_TRY_CATCH_RTN(;,STO_RTN_ERR)

	return STO_RTN_OK;
}

/** termination test */
bool DecDdMasterDSB::terminate()
{
	DSPdebugMessage("Optimality gap %e Linearization error %e |g|^2 %e\n",
			upperbound_ - valueAtProx_, phi_t_ - tau_ * gg_, gg_);

	if (upperbound_ - valueAtProx_ < eps_opt_)
		return true;

	if (phi_t_ - tau_ * gg_ < eps_e_ && gg_ < eps_g_)
		return true;

	return false;
}
