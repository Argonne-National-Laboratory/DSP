/*
 * DwBundleDual.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: kibaekkim
 */

#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#include "Utility/DspUtility.h"
#include <DantzigWolfe/DwBundleDual.h>

DwBundleDual::DwBundleDual(DwWorker* worker):
DwMaster(worker),
v_(0.0),
counter_(1),
u_(1.0),
eps_(COIN_DBL_MAX),
absp_(COIN_DBL_MAX),
alpha_(COIN_DBL_MAX),
linerr_(COIN_DBL_MAX) {
}

DwBundleDual::~DwBundleDual() {
}

DSP_RTN_CODE DwBundleDual::solve() {
	BGN_TRY_CATCH

	itercnt_ = 0;
	t_total_ = CoinGetTimeOfDay();
	t_master_ = 0.0;
	t_colgen_ = 0.0;

	/** update quadratic term */
	u_ = 1.0;
	updateCenter(u_);

	/** generate initial columns */
	DSP_RTN_CHECK_RTN_CODE(initialColumns());
	bestdualobj_ = std::min(bestdualobj_, dualobj_);

	/**
	 * The codes below are experimental to see if deactivating some dual variables would help convergence.
	 */
#if 0
	/** deactivate some dual variables (by fixed to zeros) */
	std::vector<pairIntDbl> weight;
	const CoinPackedMatrix* mat = si_->getMatrixByCol();
	for (int j = nrows_conv_; j < si_->getNumCols(); ++j) {
		const CoinShallowPackedVector col = mat->getVector(j);
		double val = 0.0;
		for (int i = 0; i < col.getNumElements(); ++i)
			val -= col.getElements()[i];
		weight.push_back(std::make_pair(j,val));
	}
	std::sort(weight.begin(), weight.end(), compPair);

	/** FIXME: Let's try 90% activation */
	int ndeactive = nrows_orig_ - floor(nrows_orig_*0.9);
	for (int j = 0; j < ndeactive; ++j) {
		//printf("Fixed column(%d) bounds to zeros.\n", weight[j].first);
		si_->setColBounds(weight[j].first, 0.0, 0.0);
	}
#endif

	DSP_RTN_CHECK_RTN_CODE(gutsOfSolve());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::createProblem() {
	BGN_TRY_CATCH

	clbd_node_ = new double [ncols_orig_];
	cubd_node_ = new double [ncols_orig_];
	CoinCopyN(clbd_orig_, ncols_orig_, clbd_node_);
	CoinCopyN(cubd_orig_, ncols_orig_, cubd_node_);

	DSP_RTN_CHECK_RTN_CODE(createPrimalProblem());
	DSP_RTN_CHECK_RTN_CODE(createDualProblem());

	/** always phase2 */
	phase_ = 2;

	/** maximization in dual */
	bestdualobj_ = COIN_DBL_MAX;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::createPrimalProblem() {
	BGN_TRY_CATCH

	/** create column-wise matrix and set number of rows */
	std::shared_ptr<CoinPackedMatrix> mat(new CoinPackedMatrix(true, 0, 0));
	mat->setDimensions(nrows_, 0);

	/** row bounds */
	std::vector<double> rlbd(nrows_, 1.0);
	std::vector<double> rubd(nrows_, 1.0);
	CoinCopyN(rlbd_orig_, nrows_orig_, &rlbd[nrows_conv_]);
	CoinCopyN(rubd_orig_, nrows_orig_, &rubd[nrows_conv_]);

	/** create solver */
	primal_si_.reset(new OsiCpxSolverInterface());

	/** load problem data */
	primal_si_->loadProblem(*mat, NULL, NULL, NULL, &rlbd[0], &rubd[0]);

	/** set display */
	primal_si_->messageHandler()->setLogLevel(0);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::createDualProblem() {
	BGN_TRY_CATCH

	/** master problem */
	std::shared_ptr<CoinPackedMatrix> mat(nullptr);
	std::vector<double> clbd(nrows_), cubd(nrows_), obj(nrows_);

	/** necessary dual variables */
	bestdualsol_.resize(nrows_);
	dualsol_.resize(nrows_);
	std::fill(bestdualsol_.begin(), bestdualsol_.begin() + nrows_conv_, COIN_DBL_MAX);
	std::fill(bestdualsol_.begin() + nrows_conv_, bestdualsol_.end(), 0.0);
	std::fill(dualsol_.begin(), dualsol_.begin() + nrows_conv_, COIN_DBL_MAX);
	std::fill(dualsol_.begin() + nrows_conv_, dualsol_.end(), 0.0);

	/** other initialization */
	p_.resize(nrows_orig_);

	/** create row-wise matrix and set number of rows */
	mat.reset(new CoinPackedMatrix(false, 0, 0));
	mat->setDimensions(0, nrows_);

	std::fill(clbd.begin(), clbd.begin() + nrows_conv_, -COIN_DBL_MAX);
	std::fill(cubd.begin(), cubd.begin() + nrows_conv_, +COIN_DBL_MAX);
	std::fill(obj.begin(), obj.begin() + nrows_conv_, -1.0);
	for (int i = 0; i < nrows_orig_; ++i) {
		clbd[nrows_conv_+i] = 0.0;
		cubd[nrows_conv_+i] = 0.0;
		if (rlbd_orig_[i] > -1.0e+20)
			cubd[nrows_conv_+i] = COIN_DBL_MAX;
		if (rubd_orig_[i] < 1.0e+20)
			clbd[nrows_conv_+i] = -COIN_DBL_MAX;
		obj[nrows_conv_+i] = -u_*bestdualsol_[nrows_conv_+i];
	}

	/** create solver */
	si_ = new OsiCpxSolverInterface();

	OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
	if (cpx) {
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_THREADS, par_->getIntParam("NUM_CORES"));
	}

	/** display */
	si_->messageHandler()->setLogLevel(std::max(-1,par_->getIntParam("LOG_LEVEL")-4));

	/** load problem data */
	si_->loadProblem(*mat, &clbd[0], &cubd[0], &obj[0], NULL, NULL);

	/** set quadratic objective term */
	updateCenter(u_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::updateCenter(double penalty) {
	OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
	if (cpx) {
		for (int j = nrows_conv_; j < nrows_; ++j) {
			si_->setObjCoeff(j, -u_*bestdualsol_[j]);
			CPXchgqpcoef(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), j, j, penalty);
		}
	} else {
		CoinError("This supports OsiCpxSolverInterface only.", "updateQuadratic", "DwBundle");
		return DSP_RTN_ERR;
	}
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::solveMaster() {
	BGN_TRY_CATCH

	OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_);
	if (!cpx) {
		CoinError("Failed to case Osi to OsiCpx", "solveMaster", "DwBundle");
		return DSP_RTN_ERR;
	}

	double stime = CoinGetTimeOfDay();
	CPXqpopt(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_PROBLEM));
	t_master_ += CoinGetTimeOfDay() - stime;
	int cpxstat = CPXgetstat(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
	switch (cpxstat) {
	case CPX_STAT_OPTIMAL:
	case CPX_STAT_NUM_BEST:
		status_ = DSP_STAT_OPTIMAL;
		break;
	case CPX_STAT_INFEASIBLE:
		status_ = DSP_STAT_PRIM_INFEASIBLE;
		break;
	case CPX_STAT_UNBOUNDED:
	case CPX_STAT_INForUNBD:
		status_ = DSP_STAT_DUAL_INFEASIBLE;
		break;
	case CPX_STAT_ABORT_OBJ_LIM:
	case CPX_STAT_ABORT_PRIM_OBJ_LIM:
		status_ = DSP_STAT_LIM_PRIM_OBJ;
		break;
	case CPX_STAT_ABORT_DUAL_OBJ_LIM:
		status_ = DSP_STAT_LIM_DUAL_OBJ;
		break;
	case CPX_STAT_ABORT_IT_LIM:
	case CPX_STAT_ABORT_TIME_LIM:
		status_ = DSP_STAT_LIM_ITERorTIME;
		break;
	case CPX_STAT_ABORT_USER:
		status_ = DSP_STAT_ABORT;
		break;
	default:
		message_->print(1, "Unexpected CPLEX status %d\n", cpxstat);
		status_ = DSP_STAT_UNKNOWN;
		break;
	}

	switch(status_) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {

		dualsol_.assign(si_->getColSolution(), si_->getColSolution() + si_->getNumCols());

		std::vector<double> d(nrows_orig_, 0.0);
		dualobj_ = si_->getObjValue();
		absp_ = 0.0;
		for (int j = 0; j < nrows_orig_; ++j) {
			d[j] = dualsol_[nrows_conv_+j] - bestdualsol_[nrows_conv_+j];
			p_[j] = -u_ * d[j];
			dualobj_ += 0.5 * u_ * bestdualsol_[nrows_conv_+j]*bestdualsol_[nrows_conv_+j];
			dualobj_ -= 0.5 * u_ * d[j] * d[j];
			absp_ += fabs(p_[j]);
		}

		v_ = dualobj_ - bestdualobj_;

		/** adjust v if subproblem was not solved to optimality */
		for (auto it = status_subs_.begin(); it != status_subs_.end(); it++)
			if (*it != DSP_STAT_OPTIMAL) {
				v_ = std::min(v_, -1.0e-4);
				break;
			}

		alpha_ = -v_;
		for (int j = 0; j < nrows_orig_; ++j)
			alpha_ -= p_[j] * d[j];

		/** get primal solutions by solving the primal master */
		primal_si_->resolve();
		if (primal_si_->isProvenOptimal()) {
			primobj_ = primal_si_->getObjValue();
			primsol_.assign(primal_si_->getColSolution(), primal_si_->getColSolution() + primal_si_->getNumCols());
		} else {
			message_->print(3, "  The primal master could not be solved to optimality.\n");
			primobj_ = COIN_DBL_MAX;
			primsol_.assign(si_->getRowPrice(), si_->getRowPrice() + si_->getNumRows());
		}

		absgap_ = fabs(primobj_+bestdualobj_);
		relgap_ = absgap_/(1.0e-10+fabs(bestdualobj_));
		break;
	}
	default:
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::updateModel() {
	double u = u_, newu;

	BGN_TRY_CATCH

	/** descent test */
	bool foundBetter = false;
	if (dualobj_ <= bestdualobj_ + mL_ * v_) {
		message_->print(2, "  Serious step: best dual %+e -> %e\n", -bestdualobj_, -dualobj_);
		foundBetter = true;
		/** reset subproblem time increment */
		worker_->resetTimeIncrement();
	}

	/** update weight */
	if (foundBetter) {
		if (dualobj_ <= bestdualobj_ + mR_ * v_ && counter_ > 0)
			u = 2 * u_ * (1 - (dualobj_ - bestdualobj_) / v_);
		else if (counter_ > 3)
			u = 0.5 * u_;
		newu = std::max(std::max(u, 0.1*u_), umin_);
		eps_ = std::max(eps_, -2*v_);
		counter_ = std::max(counter_+1,1);
		if (u_ != newu)
			counter_ = 1;
		u_ = newu;
		bestdualobj_ = dualobj_;
		bestdualsol_ = dualsol_;
	} else {
		eps_ = std::min(eps_, absp_ + alpha_);
		if (linerr_ > std::max(eps_, -10*v_) && counter_ < -3)
			u = 2 * u_ * (1 - (dualobj_ - bestdualobj_) / v_);
		newu = std::min(u, 10*u_);
		counter_ = std::min(counter_-1,-1);
		if (u_ != newu)
			counter_ = -1;
		u_ = newu;
	}
	updateCenter(u_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwBundleDual::terminationTest(int nnewcols) {
	BGN_TRY_CATCH

	if (relgap_ <= 1.0e-4)
		return true;

	if (v_ >= 0.0)
		return true;

	if (iterlim_ <= itercnt_)
		return true;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return false;
}

DSP_RTN_CODE DwBundleDual::addCols(
		const double* piA,                    /**< [in] pi^T A */
		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
		std::vector<int>& statuses,           /**< [in] subproblem solution status */
		std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
		std::vector<double>& objs,            /**< [in] subproblem objective values */
		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
		int& nadded                           /**< [out] number of columns added */) {
	return addRows(indices, statuses, cxs, objs, sols, nadded);
}

DSP_RTN_CODE DwBundleDual::addRows(
		std::vector<int>& indices,            /**< [in] subproblem indices corresponding to cols*/
		std::vector<int>& statuses,           /**< [in] subproblem solution status */
		std::vector<double>& cxs,             /**< [in] solution times original objective coefficients */
		std::vector<double>& objs,            /**< [in] subproblem objective values */
		std::vector<CoinPackedVector*>& sols, /**< [in] subproblem solutions */
		int& nadded                           /**< [out] number of columns added */) {
	CoinPackedVector cutvec;
	double cutcoef, cutrhs;

	BGN_TRY_CATCH

	/** allocate memory */
	std::vector<double> Ax(nrows_orig_, 0.0);
	cutvec.reserve(nrows_conv_+nrows_orig_);

	/** reset counter */
	nadded = 0;

	linerr_ = bestdualobj_;

	for (unsigned int s = 0; s < indices.size(); ++s) {
		int sind = indices[s]; /**< actual subproblem index */

		/** retrieve subproblem solution */
		const CoinPackedVector* x = sols[s];
#if 0
		for (int i = 0; i < x->getNumElements(); ++i) {
			int j = x->getIndices()[i];
			if (j >= ncols_orig_) break;
			if (x->getElements()[i] > cubd_node_[j] ||
				x->getElements()[i] < clbd_node_[j])
				printf("violation in subproblem %d (col %d, status %d): %+e <= %+e <= %+e\n",
						sind, j, statuses[s], clbd_node_[j], x->getElements()[i], cubd_node_[j]);
		}
#endif

		/** take A x^k */
		mat_orig_->times(*x, &Ax[0]);

		/** dual objective and linearization error */
		linerr_ -= cxs[s];
		for (int i = 0; i < nrows_orig_; ++i) {
			linerr_ += Ax[i] * bestdualsol_[nrows_conv_+i];
			if (rlbd_orig_[i] > -1.0e+20)
				linerr_ -= bestdualsol_[nrows_conv_+i] * rlbd_orig_[i];
			if (rubd_orig_[i] < 1.0e+20)
				linerr_ -= bestdualsol_[nrows_conv_+i] * rubd_orig_[i];
		}

		if (dualsol_[sind] > objs[s] + 1.0e-6 ||
				statuses[s] == DSP_STAT_DUAL_INFEASIBLE) {

			/** clear cut vector */
			cutvec.clear();

			/** original constraints */
			if (statuses[s] != DSP_STAT_DUAL_INFEASIBLE)
				cutvec.insert(sind, 1.0);
			cutrhs = cxs[s];
			for (int i = 0; i < nrows_orig_; ++i) {
				cutcoef = Ax[i];
				if (rlbd_orig_[i] > -1.0e+20)
					cutcoef -= rlbd_orig_[i];
				if (rubd_orig_[i] < 1.0e+20)
					cutcoef -= rubd_orig_[i];
				if (fabs(cutcoef) > 1.0e-10)
					cutvec.insert(nrows_conv_+i, cutcoef);
			}

			/** add row to the dual master */
			si_->addRow(cutvec, -COIN_DBL_MAX, cutrhs);

			/** add column to the primal master */
			primal_si_->addCol(cutvec, 0.0, COIN_DBL_MAX, cutrhs);

			/** store columns */
			cols_generated_.push_back(new DwCol(sind, *x, cutvec, cutrhs, 0.0, COIN_DBL_MAX));
			nadded++;
		}
	}
	DSPdebugMessage("Number of columns in the pool: %u\n", cols_generated_.size());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::getLagrangianBound(
		std::vector<double>& objs /**< [in] subproblem objective values */) {
	BGN_TRY_CATCH

	/** calculate lower bound */
	dualobj_ = 0.0;
	for (int j = 0; j < nrows_orig_; ++j) {
		if (rlbd_orig_[j] > -1.0e+20)
			dualobj_ -= dualsol_[nrows_conv_+j] * rlbd_orig_[j];
		else if (rubd_orig_[j] < 1.0e+20)
			dualobj_ -= dualsol_[nrows_conv_+j] * rubd_orig_[j];
	}
	for (unsigned int s = 0; s < objs.size(); ++s)
		dualobj_ -= objs[s];

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

void DwBundleDual::setBranchingObjects(const DspBranch* branchobj) {
	/** shouldn't be null */
	if (branchobj == NULL)
		return;

	BGN_TRY_CATCH

	/** restore column bounds */
	CoinCopyN(clbd_orig_, ncols_orig_, clbd_node_);
	CoinCopyN(cubd_orig_, ncols_orig_, cubd_node_);

	/** update column bounds at the current node */
	for (unsigned j = 0; j < branchobj->index_.size(); ++j) {
		clbd_node_[branchobj->index_[j]] = branchobj->lb_[j];
		cubd_node_[branchobj->index_[j]] = branchobj->ub_[j];
		message_->print(3, "  Branching index %d: [%+e, %+e]\n", branchobj->index_[j], branchobj->lb_[j], branchobj->ub_[j]);
	}

	/** remove all columns in the primal master */
	std::vector<int> delcols(primal_si_->getNumCols());
	CoinIotaN(&delcols[0], primal_si_->getNumCols(), 0);
	primal_si_->deleteCols(primal_si_->getNumCols(), &delcols[0]);

	/** remove all rows in the dual master */
	std::vector<int> delrows(si_->getNumRows());
	CoinIotaN(&delrows[0], si_->getNumRows(), 0);
	si_->deleteRows(si_->getNumRows(), &delrows[0]);
#if 1
	/** add branching rows */
	for (auto it = cols_generated_.begin(); it != cols_generated_.end(); it++) {
		(*it)->active_ = true;
		for (unsigned j = 0, i = 0; j < branchobj->index_.size(); ++j) {
			int sparse_index = (*it)->x_.findIndex(branchobj->index_[j]);
			double val = 0.0;
			if (sparse_index > -1)
				val = (*it)->x_.getElements()[sparse_index];
			if (val < branchobj->lb_[j] || val > branchobj->ub_[j]) {
				(*it)->active_ = false;
				break;
			}
		}

		if ((*it)->active_) {
			si_->addRow((*it)->col_, -COIN_DBL_MAX, (*it)->obj_);
			primal_si_->addCol((*it)->col_, 0.0, COIN_DBL_MAX, (*it)->obj_);
		}
	}
#else
	for (auto it = cols_generated_.begin(); it != cols_generated_.end(); it++)
		(*it)->active_ = false;
#endif
	message_->print(3, "Appended dynamic columns in the master (%d / %u cols).\n", si_->getNumRows(), cols_generated_.size());

	/** apply column bounds */
	std::vector<int> ncols_inds(ncols_orig_);
	CoinIotaN(&ncols_inds[0], ncols_orig_, 0);
	worker_->setColBounds(ncols_orig_, &ncols_inds[0], clbd_node_, cubd_node_);

	/** set known best bound */
	bestdualobj_ = branchobj->bestBound_;

	END_TRY_CATCH(;)
}

void DwBundleDual::printIterInfo() {
	message_->print(1, "Iteration %3d: Primal %+e, Best %+e ",
			itercnt_, primobj_, -bestdualobj_);
	if (relgap_ < 1000)
		message_->print(1, "(gap %.2f %%), ", relgap_*100);
	else
		message_->print(1, "(gap Large %%), ");
	message_->print(1, "nrows %d, ncols %d, ", si_->getNumRows(), si_->getNumCols());
	message_->print(1, "timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
			CoinGetTimeOfDay() - t_total_, t_master_, t_colgen_, status_);
	message_->print(3, "  predicted ascent %+e, |p| %+e, alpha %+e, linerr %+e, eps %+e, u %+e, counter %d\n",
			-v_, absp_, alpha_, -linerr_, eps_, u_, counter_);
}
