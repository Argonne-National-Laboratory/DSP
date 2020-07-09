/*
 * DwBundleDual.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include "SolverInterface/DspOsiScip.h"
#include "SolverInterface/DspOsiCpx.h"
#include "SolverInterface/DspOsiGrb.h"
#include "Solver/DantzigWolfe/DwBundleDual.h"
#include "Utility/DspUtility.h"

DwBundleDual::DwBundleDual(DwWorker* worker):
DwMaster(worker),
v_(0.0),
counter_(1),
u_(1.0),
eps_(COIN_DBL_MAX),
absp_(COIN_DBL_MAX),
alpha_(COIN_DBL_MAX),
linerr_(COIN_DBL_MAX),
prev_dualobj_(COIN_DBL_MAX),
nstalls_(0),
numFixedRows_(0) {
}

DwBundleDual::DwBundleDual(const DwBundleDual& rhs):
DwMaster(rhs),
v_(rhs.v_),
counter_(rhs.counter_),
u_(rhs.u_),
eps_(rhs.eps_),
absp_(rhs.absp_),
alpha_(rhs.alpha_),
linerr_(rhs.linerr_),
prev_dualobj_(rhs.prev_dualobj_),
nstalls_(rhs.nstalls_),
numFixedRows_(rhs.numFixedRows_) {
	d_ = rhs.d_;
	p_ = rhs.p_;
	primal_si_.reset(rhs.primal_si_->clone());
}

DwBundleDual& DwBundleDual::operator =(const DwBundleDual& rhs) {
	DwMaster::operator =(rhs);
	v_ = rhs.v_;
	counter_ = rhs.counter_;
	u_ = rhs.u_;
	eps_ = rhs.eps_;
	absp_ = rhs.absp_;
	alpha_ = rhs.alpha_;
	linerr_ = rhs.linerr_;
	prev_dualobj_ = rhs.prev_dualobj_;
	nstalls_ = rhs.nstalls_;
	numFixedRows_ = rhs.numFixedRows_;
	d_ = rhs.d_;
	p_ = rhs.p_;
	primal_si_.reset(rhs.primal_si_->clone());
	return *this;
}

DwBundleDual::~DwBundleDual() {
}

DSP_RTN_CODE DwBundleDual::solve() {
	BGN_TRY_CATCH

	itercnt_ = 0;
	t_start_ = CoinGetTimeOfDay();
	t_total_ = 0.0;
	t_master_ = 0.0;
	t_colgen_ = 0.0;
	status_ = DSP_STAT_FEASIBLE;

	/** clear logs */
	log_time_.clear();
	log_bestdual_bounds_.clear();

	/** initial price to generate columns */
	bestdualsol_.assign(nrows_, 0.0);
	dualsol_ = bestdualsol_;
	std::fill(dualsol_.begin(), dualsol_.begin() + nrows_conv_, COIN_DBL_MAX);

	/** update quadratic term */
	u_ = par_->getDblParam("DW/INIT_CENTER");
	counter_ = 0;
	eps_ = COIN_DBL_MAX;
	DSP_RTN_CHECK_THROW(updateCenter(u_));

	DSPdebugMessage("Initial Setting at solve()\n");

	if (osi_ == NULL || getSiPtr()->getNumRows() == numFixedRows_) {
		/** generate initial columns */
		DSPdebugMessage("Generating columns..\n");
		double stime = CoinGetTimeOfDay();
		DSP_RTN_CHECK_THROW(generateCols());
		t_colgen_ += CoinGetTimeOfDay() - stime;
		DSPdebugMessage("Generated columns..\n");

		/** subproblem solution may declare infeasibility. */
		for (auto st = status_subs_.begin(); st != status_subs_.end(); st++)
			if (*st == DSP_STAT_PRIM_INFEASIBLE) {
				status_ = DSP_STAT_PRIM_INFEASIBLE;
				message_->print(1, "Subproblem solution is infeasible.\n");
				break;
			}

		/** check time limit */
		t_total_ = CoinGetTimeOfDay() - t_start_;
		if (time_remains_ < t_total_) {
			message_->print(1, "Time limit reached.\n");
			status_ = DSP_STAT_LIM_ITERorTIME;
		}

		if (status_ == DSP_STAT_FEASIBLE) {
			message_->print(1, "Generated %u initial columns. Initial dual bound %.12e\n", ngenerated_, -dualobj_);
			if (dualobj_ < bestdualobj_) {
				bestdualobj_ = dualobj_;
				bestdualsol_ = dualsol_;
			}
			DSP_RTN_CHECK_THROW(gutsOfSolve());
		}		
	} else {
		message_->print(1, "Starting with %u existing columns.\n", getSiPtr()->getNumRows());
		DSP_RTN_CHECK_THROW(gutsOfSolve());
	}


	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::createProblem() {
	BGN_TRY_CATCH

	clbd_node_ = clbd_orig_;
	cubd_node_ = cubd_orig_;

	DSP_RTN_CHECK_THROW(createPrimalProblem());
	DSP_RTN_CHECK_THROW(createDualProblem());

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
	std::vector<double> rlbd(nrows_);
	std::vector<double> rubd(nrows_);
	std::fill(rlbd.begin(), rlbd.begin() + nrows_conv_, 1.0);
	std::fill(rubd.begin(), rubd.begin() + nrows_conv_, 1.0);
	std::copy(rlbd_orig_.begin(), rlbd_orig_.begin() + nrows_orig_, rlbd.begin() + nrows_conv_);
	std::copy(rubd_orig_.begin(), rubd_orig_.begin() + nrows_orig_, rubd.begin() + nrows_conv_);

	/** create solver */
	DspOsi * osi = createDspOsi();
	if (osi)
		primal_si_.reset(osi);
	else
		throw CoinError("Failed to create DspOsi", "createPrimalProblem", "DwBundleDual");

	/** load problem data */
	primal_si_->si_->loadProblem(*mat, NULL, NULL, NULL, &rlbd[0], &rubd[0]);

	/** set display */
	primal_si_->setLogLevel(par_->getIntParam("DW/MASTER/SOLVER/LOG_LEVEL"));

	/** set number of cores */
	primal_si_->setNumCores(par_->getIntParam("NUM_CORES"));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

void DwBundleDual::initDualSolver(
		const CoinPackedMatrix& m, 
		std::vector<double>& clbd, 
		std::vector<double>& cubd, 
		std::vector<double>& obj, 
		std::vector<double>& rlbd, 
		std::vector<double>& rubd) {
	BGN_TRY_CATCH

	/** create solver */
	osi_ = createDspOsi();

	/** set display */
	osi_->setLogLevel(par_->getIntParam("DW/MASTER/SOLVER/LOG_LEVEL"));

	/** load problem data */
	if (rlbd.size() == 0)
		getSiPtr()->loadProblem(m, &clbd[0], &cubd[0], &obj[0], NULL, NULL);
	else
		getSiPtr()->loadProblem(m, &clbd[0], &cubd[0], &obj[0], &rlbd[0], &rubd[0]);

	END_TRY_CATCH(;)
}

DSP_RTN_CODE DwBundleDual::createDualProblem() {
	BGN_TRY_CATCH

	/** master problem */
	std::shared_ptr<CoinPackedMatrix> mat(nullptr);
	std::vector<double> clbd(nrows_), cubd(nrows_), obj(nrows_);

	/** necessary dual variables */
	dualsol_.resize(nrows_);
	std::fill(dualsol_.begin(), dualsol_.begin() + nrows_conv_, COIN_DBL_MAX);
	std::fill(dualsol_.begin() + nrows_conv_, dualsol_.end(), 0.0);
	bestdualsol_ = dualsol_;

	/** other initialization */
	p_.reserve(nrows_orig_+nrows_branch_);
	d_.reserve(nrows_orig_+nrows_branch_);
	std::fill(p_.begin(), p_.end(), 0.0);
	std::fill(d_.begin(), d_.end(), 0.0);

	/** create row-wise matrix and set number of rows */
	mat.reset(new CoinPackedMatrix(false, 0, 0));
	mat->setDimensions(0, nrows_);

	std::fill(clbd.begin(), clbd.begin() + nrows_conv_, -COIN_DBL_MAX);
	std::fill(cubd.begin(), cubd.begin() + nrows_conv_, +COIN_DBL_MAX);
	std::fill(obj.begin(), obj.begin() + nrows_conv_, -1.0);
	for (int i = 0; i < nrows_orig_; ++i) {
		clbd[nrows_conv_+i] = 0.0;
		cubd[nrows_conv_+i] = 0.0;
		obj[nrows_conv_+i] = -u_*bestdualsol_[nrows_conv_+i];
		if (rlbd_orig_[i] > -1.0e+20) {
			cubd[nrows_conv_+i] = COIN_DBL_MAX;
			obj[nrows_conv_+i] -= rlbd_orig_[i];
		}
		if (rubd_orig_[i] < 1.0e+20) {
			clbd[nrows_conv_+i] = -COIN_DBL_MAX;
			obj[nrows_conv_+i] -= rubd_orig_[i];
		}
	}

	/** initialize external solver */
	std::vector<double> rlbd, rubd;
	initDualSolver(*mat, clbd, cubd, obj, rlbd, rubd);

	/** set quadratic objective term */
	DSP_RTN_CHECK_THROW(updateCenter(u_));

	/** set number of cores */
	osi_->setNumCores(par_->getIntParam("NUM_CORES"));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::updateCenter(double penalty) {
#define FREE_MEMORY \
	FREE_PTR(qmat)

	CoinPackedMatrix *qmat = NULL;

	BGN_TRY_CATCH

	double coef;
	const double* rlbd = primal_si_->si_->getRowLower();
	const double* rubd = primal_si_->si_->getRowUpper();

	/** update linear coefficients */
	for (int j = nrows_conv_; j < nrows_; ++j) {
		coef = -penalty*bestdualsol_[j];
		if (fabs(rlbd[j]) < 1.0e+20)
			coef -= rlbd[j];
		if (fabs(rubd[j]) < 1.0e+20)
			coef -= rubd[j];
		if (fabs(coef) >= 1.0e+20) {
			printf("penalty = %e, bestdualsol_[%d] = %e, coef = %e\n", penalty, j, bestdualsol_[j], coef);
			CoinError("Invalid objective coefficient.", "updateCenter", "DwBundleDual");
			return DSP_RTN_ERR;
		}
		getSiPtr()->setObjCoeff(j, coef);
	}

	/** update quadratic term */
	vector<int> indices;
	indices.reserve(nrows_ - nrows_conv_);
	for (int j = nrows_conv_; j < nrows_; ++j) {
		indices.push_back(j);
	}
	vector<double> elements(nrows_ - nrows_conv_, penalty);

	qmat = new CoinPackedMatrix(true, &indices[0], &indices[0], &elements[0], nrows_-nrows_conv_);
	osi_->loadQuadraticObjective(*qmat);

	FREE_MEMORY

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::callMasterSolver() {

	if (par_->getIntParam("DW/MASTER/SOLVER") == OsiCpx) {
#ifdef DSP_HAS_CPX

		OsiCpxSolverInterface* cpx = NULL;

		BGN_TRY_CATCH

		cpx = dynamic_cast<OsiCpxSolverInterface*>(getSiPtr());

		END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

		/** use dual simplex for QP, which is numerically much more stable than Barrier */
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_QPMETHOD, 0);
		CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_NUMERICALEMPHASIS, 0);

		int pass = 0;
		while (pass <= 2) {
			CPXqpopt(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_PROBLEM));
			int cpxstat = CPXgetstat(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
			message_->print(5, "CPLEX status %d\n", cpxstat);
			switch (cpxstat) {
			case CPX_STAT_OPTIMAL:
			case CPX_STAT_NUM_BEST:
			case CPX_STAT_OPTIMAL_INFEAS:
				status_ = DSP_STAT_OPTIMAL;
				pass = 10;
				break;
			case CPX_STAT_INFEASIBLE:
				status_ = DSP_STAT_DUAL_INFEASIBLE;
				message_->print(0, "Unexpected CPLEX status %d\n", cpxstat);
				pass = 10;
				break;
			case CPX_STAT_UNBOUNDED:
			case CPX_STAT_INForUNBD:
				if (pass == 0) {
					CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_NUMERICALEMPHASIS, 1);
					message_->print(1, "Re-running CPLEX with numerical emphasis\n");
				} else if (pass == 1) {
					CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_QPMETHOD, 2);
					message_->print(1, "Re-running CPLEX with dual simplex\n");
				} else {
					status_ = DSP_STAT_PRIM_INFEASIBLE;
					message_->print(1, "CPLEX failed to optimize\n");
					message_->print(0, "Unexpected CPLEX status %d\n", cpxstat);
					getSiPtr()->writeMps("unbounded_master");
				}
				pass++;
				break;
			case CPX_STAT_ABORT_OBJ_LIM:
			case CPX_STAT_ABORT_PRIM_OBJ_LIM:
				if (pass == 0) {
					CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_NUMERICALEMPHASIS, 1);
					message_->print(1, "Re-running CPLEX with numerical emphasis\n");
				} else if (pass == 1) {
					CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_QPMETHOD, 2);
					message_->print(1, "Re-running CPLEX with dual simplex\n");
				} else {
					status_ = DSP_STAT_LIM_PRIM_OBJ;
					message_->print(1, "CPLEX failed to optimize\n");
					message_->print(0, "Unexpected CPLEX status %d\n", cpxstat);
					getSiPtr()->writeMps("unbounded_master");
				}
				pass++;
				break;
			case CPX_STAT_ABORT_DUAL_OBJ_LIM:
				status_ = DSP_STAT_LIM_DUAL_OBJ;
				pass = 10;
				break;
			case CPX_STAT_ABORT_IT_LIM:
			case CPX_STAT_ABORT_TIME_LIM:
				status_ = DSP_STAT_LIM_ITERorTIME;
				pass = 10;
				break;
			case CPX_STAT_ABORT_USER:
				status_ = DSP_STAT_ABORT;
				pass = 10;
				break;
			default:
				message_->print(0, "Unexpected CPLEX status %d\n", cpxstat);
				status_ = DSP_STAT_UNKNOWN;
				pass = 10;
				break;
			}
		}
#endif
	} else {
		osi_->solve();
		status_ = osi_->status();
	}
	return status_;
}

void DwBundleDual::assignMasterSolution(std::vector<double>& sol) {
	sol.assign(getSiPtr()->getColSolution(), getSiPtr()->getColSolution() + getSiPtr()->getNumCols());
}

DSP_RTN_CODE DwBundleDual::solveMaster() {
	BGN_TRY_CATCH

	/** get previous objective */
	prev_dualobj_ = dualobj_;

	/** call solver */
	status_ = callMasterSolver();

	/** status_ must be set at this point. */
	DSPdebugMessage("master status_ %d\n", status_);

	switch(status_) {
	case DSP_STAT_OPTIMAL:
	case DSP_STAT_FEASIBLE:
	case DSP_STAT_LIM_ITERorTIME: {

		status_ = DSP_STAT_OPTIMAL;

		assignMasterSolution(dualsol_);
#ifdef DSP_DEBUG
		printf("bestdualsol_:\n");
		DspMessage::printArray(bestdualsol_.size(), &bestdualsol_[0]);
		printf("dualsol_:\n");
		DspMessage::printArray(dualsol_.size(), &dualsol_[0]);
#endif

		d_.resize(nrows_-nrows_conv_);
		p_.resize(nrows_-nrows_conv_);
		double polyapprox = 0.0;
		absp_ = 0.0;
		for (int j = 0; j < nrows_conv_; ++j)
			polyapprox -= getSiPtr()->getColSolution()[j];
		for (int j = nrows_conv_; j < nrows_; ++j) {
			d_[j-nrows_conv_] = dualsol_[j] - bestdualsol_[j];
			p_[j-nrows_conv_] = -u_ * d_[j-nrows_conv_];
			// polyapprox += u_ * dualsol_[j] * (bestdualsol_[j] - 0.5 * dualsol_[j]);
			absp_ += fabs(p_[j-nrows_conv_]);
			//printf("j %d: dualsol %+e, bestdualsol %+e, d %+e, p %+e, polyapprox %+e\n",
			//		j, dualsol_[j], bestdualsol_[j], d_[j-nrows_conv_], p_[j-nrows_conv_], polyapprox);
		}
		//printf("u*dualsol*(bestdualsol-0.5*dualsol) = %+e\n", polyapprox);
		// polyapprox += getObjValue();

		DSPdebugMessage("getObjValue %e, polyapprox %e, bestdualobj_ %e\n", getObjValue(), polyapprox, bestdualobj_);
		v_ = polyapprox - bestdualobj_;

		/** adjust v if subproblem was not solved to optimality */
		for (auto it = status_subs_.begin(); it != status_subs_.end(); it++)
			if (*it != DSP_STAT_OPTIMAL) {
				v_ = std::min(v_, -1.0e-4);
				break;
			}

		alpha_ = -v_;
		for (int j = nrows_conv_; j < nrows_; ++j)
			alpha_ += p_[j-nrows_conv_] * d_[j-nrows_conv_];

		/** get primal solutions by solving the primal master */
		DSPdebug(primal_si_->si_->writeMps("PrimMaster"));
		primal_si_->si_->resolve();
		if (primal_si_->si_->isProvenOptimal()) {
			primobj_ = primal_si_->getPrimObjValue();
			primsol_.assign(primal_si_->si_->getColSolution(), primal_si_->si_->getColSolution() + primal_si_->si_->getNumCols());
		} else {
			message_->print(5, "  The primal master could not be solved to optimality.\n");
			primobj_ = COIN_DBL_MAX;
			//primsol_.assign(si_->getRowPrice(), si_->getRowPrice() + si_->getNumRows());
		}

		break;
	}
	default:
		char msg[128];
		sprintf(msg, "Unexpected solution status: %d", status_);
		getSiPtr()->writeMps("unexpected_master");
		throw CoinError(msg, "solveMaster", "DwBundleDual");
		break;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::updateModel() {
	BGN_TRY_CATCH

	double u = u_, newu;

	/** descent test */
	bool foundBetter = false;
	DSPdebugMessage("dualobj_ %e bestdualobj_ %e serious step? %e\n", dualobj_, bestdualobj_, bestdualobj_ + mL_ * v_ - dualobj_);
	if (dualobj_ <= bestdualobj_ + mL_ * v_ && v_ < 0) {
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
		newu = std::min(std::max(std::max(u, 0.1*u_), umin_), umax_);
		eps_ = std::max(eps_, -2*v_);

		/** update counter */
		if (fabs(u_-newu) > 1.0e-8)
			counter_ = 1;
		else
			counter_ = std::max(counter_+1,1);

		bestdualobj_ = dualobj_;
		bestdualsol_ = dualsol_;
		nstalls_ = 0;
	} else {
		eps_ = std::min(eps_, absp_ + alpha_);
		if (-linerr_ > std::max(eps_, -10*v_) && counter_ < -3)
			u = 2 * u_ * (1 - (dualobj_ - bestdualobj_) / v_);
		newu = std::min(u, 10*u_);

		// My customization
		nstalls_ = fabs(dualobj_ - prev_dualobj_) < 1.0e-6 ? nstalls_ + 1 : 0;
		if (nstalls_ > 0)
			message_->print(3, "number of stalls: %d\n", nstalls_);
		if (nstalls_ > 3)
			newu = 0.1*u_;
		else if ((primobj_ >= 1.0e+20 && v_ >= -par_->getDblParam("DW/MIN_INCREASE")) || nstalls_ > 0)
			newu = 0.5*u_;
		// else if (counter_ < -5)
		// 	newu = 10*u_;

		/** update counter */
		if (fabs(u_-newu) > 1.0e-8)
			counter_ = -1;
		else
			counter_ = std::min(counter_-1,-1);
	}
	u_ = newu;
	DSP_RTN_CHECK_THROW(updateCenter(u_));

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

bool DwBundleDual::terminationTest() {
	BGN_TRY_CATCH

	if (-bestdualobj_ >= bestprimobj_) {
		message_->print(3, "Terminated due to best dual bound (%e) >= best primal bound (%e)\n", -bestdualobj_, bestprimobj_);
		status_ = DSP_STAT_LIM_DUAL_OBJ;
		return true;
	}

	if (primobj_ < 1.0e+20 && getRelApproxGap() <= par_->getDblParam("DW/GAPTOL"))
		return true;

	if (primobj_ < 1.0e+20 && v_ >= -par_->getDblParam("DW/MIN_INCREASE"))
		return true;

	if (iterlim_ <= itercnt_) {
		message_->print(3, "Warning: Iteration limit reached.\n");
		status_ = DSP_STAT_LIM_ITERorTIME;
		return true;
	}

	if (nstalls_ >= 30) {
		message_->print(3, "Warning: Maximum stalling limit reached.\n");
		status_ = DSP_STAT_LIM_ITERorTIME;
		return true;
	}

	if (time_remains_ < t_total_) {
		message_->print(3, "Warning: Time limit reached.\n");
		status_ = DSP_STAT_LIM_ITERorTIME;
		return true;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return false;
}

DSP_RTN_CODE DwBundleDual::addCols(
		std::vector<int>& indices,           /**< [in] subproblem indices corresponding to cols*/
		std::vector<int>& statuses,          /**< [in] subproblem solution status */
		std::vector<double>& cxs,            /**< [in] solution times original objective coefficients */
		std::vector<double>& objs,           /**< [in] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [in] subproblem solutions */) {
	return addRows(indices, statuses, cxs, objs, sols);
}

DSP_RTN_CODE DwBundleDual::addRows(
		std::vector<int>& indices,           /**< [in] subproblem indices corresponding to cols*/
		std::vector<int>& statuses,          /**< [in] subproblem solution status */
		std::vector<double>& cxs,            /**< [in] solution times original objective coefficients */
		std::vector<double>& objs,           /**< [in] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [in] subproblem solutions */) {
	CoinPackedVector cutvec;
	double cutcoef, cutrhs;

	BGN_TRY_CATCH

	/** allocate memory */
	std::vector<double> Ax(nrows_orig_, 0.0);
	cutvec.reserve(nrows_);
	assert(mat_orig_->getNumRows() == Ax.size());

	/** reset counter */
	ngenerated_ = 0;

	/** linearization error */
	linerr_ = 0.0;
	for (int i = nrows_conv_; i < nrows_; ++i) {
		if (primal_si_->si_->getRowLower()[i] > -1.0e+20)
			linerr_ += primal_si_->si_->getRowLower()[i] * d_[i-nrows_conv_];
		if (primal_si_->si_->getRowUpper()[i] < 1.0e+20)
			linerr_ += primal_si_->si_->getRowUpper()[i] * d_[i-nrows_conv_];
	}

	for (unsigned int s = 0; s < indices.size(); ++s) {
		int sind = indices[s]; /**< actual subproblem index */

		/** retrieve subproblem solution */
		const CoinPackedVector* x = sols[s];
#ifdef DSP_DEBUG
		DSPdebugMessage("sol[%d] =\n", s);
		DspMessage::printArray(x);
#endif
#if 0
		for (int i = 0; i < x->getNumElements(); ++i) {
			int j = x->getIndices()[i];
			if (j >= ncols_orig_) break;
			if (x->getElements()[i] - cubd_node_[j] > 1.0e-6 ||
				x->getElements()[i] - clbd_node_[j] < -1.0e-6)
				printf("violation in subproblem %d (col %d, status %d): %+e <= %+e <= %+e\n",
						sind, j, statuses[s], clbd_node_[j], x->getElements()[i], cubd_node_[j]);
		}
#endif

		/** take A x^k */
		mat_orig_->times(*x, &Ax[0]);

		/** clear cut vector */
		cutvec.clear();

		/** original constraints */
		if (statuses[s] != DSP_STAT_DUAL_INFEASIBLE)
			cutvec.insert(sind, 1.0);
		for (int i = 0; i < nrows_orig_; ++i) {
			cutcoef = Ax[i];
			assert(fabs(cutcoef) < 1e+20);
			linerr_ -= cutcoef * d_[i];
			if (fabs(cutcoef) > 1.0e-10)
				cutvec.insert(nrows_conv_+i, cutcoef);
		}
		cutrhs = cxs[s];

		/** branching objects */
#ifdef USE_ROW_TO_COL
		for (int i = 0; i < nrows_branch_; ++i) {
			int j = branch_row_to_col_[nrows_core_ + i];
			int sparse_index = x->findIndex(j);
			if (sparse_index == -1) continue;
			cutcoef = x->getElements()[sparse_index];
			linerr_ -= cutcoef * d_[nrows_orig_+i];
			if (fabs(cutcoef) > 1.0e-10)
				cutvec.insert(nrows_core_+i, cutcoef);
		}
#else
		if (nrows_branch_ > 0) {
			double* dx = x->denseVector(ncols_orig_);
			for (int i = 0; i < nrows_branch_; ++i) {
				cutcoef = branch_row_to_vec_[nrows_core_+i].dotProduct(dx);
				assert(fabs(cutcoef) < 1e+20);
				linerr_ -= cutcoef * d_[nrows_orig_+i];
				if (fabs(cutcoef) > 1.0e-10)
					cutvec.insert(nrows_core_+i, cutcoef);
			}
			// free dense vector
			FREE_ARRAY_PTR(dx);
		}
#endif

		DSPdebugMessage("subproblem %d: status %d, objective %+e, violation %+e\n", sind, statuses[s], objs[s], dualsol_[sind] - objs[s]);
		if (statuses[s] == DSP_STAT_UNKNOWN)
			throw "Unknown subproblem solution status";

		if (dualsol_[sind] > objs[s] + 1.0e-8 ||
				statuses[s] == DSP_STAT_DUAL_INFEASIBLE) {
#ifdef DSP_DEBUG
			DspMessage::printArray(&cutvec);
			printf("cutrhs %+e\n", cutrhs);
#endif

			/** add row to the dual master */
			addDualRow(cutvec, -COIN_DBL_MAX, cutrhs);

			/** add column to the primal master */
			primal_si_->si_->addCol(cutvec, 0.0, COIN_DBL_MAX, cutrhs);

			/** store columns */
			cols_generated_.push_back(new DwCol(sind, primal_si_->si_->getNumCols() - 1, *x, cutvec, cutrhs, 0.0, COIN_DBL_MAX));
			ngenerated_++;
		} else {
			/** store columns */
			cols_generated_.push_back(new DwCol(sind, -1, *x, cutvec, cutrhs, 0.0, COIN_DBL_MAX, false));
		}

		// free dense vector
		// FREE_ARRAY_PTR(dx);
	}
	DSPdebugMessage("Number of columns in the pool: %lu\n", cols_generated_.size());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::getLagrangianBound(
		std::vector<double>& objs /**< [in] subproblem objective values */) {
	BGN_TRY_CATCH

	/** calculate lower bound */
	dualobj_ = -std::accumulate(objs.begin(), objs.end(), 0.0);
	DSPdebugMessage("Sum of subproblem objectives: %+e\n", -dualobj_);
	const double* rlbd = primal_si_->si_->getRowLower();
	const double* rubd = primal_si_->si_->getRowUpper();
	for (int j = nrows_conv_; j < nrows_; ++j) {
		if (rlbd[j] > -1.0e+20)
			dualobj_ -= dualsol_[j] * rlbd[j];
		else if (rubd[j] < 1.0e+20)
			dualobj_ -= dualsol_[j] * rubd[j];
	}

	/** linearization error */
	linerr_ += bestdualobj_ - dualobj_;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::reduceCols(int &num_removed) {
	
	//if (nstalls_ < 3 || si_->getNumRows() <= si_->getNumCols())
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	/** row indices to delete */
	std::vector<int> delrows;
	std::vector<int> delblk(model_->getNumSubproblems(), 0);

	for (unsigned k = 0; k < cols_generated_.size(); ++k) {
		if (cols_generated_[k]->active_)
			delblk[cols_generated_[k]->blockid_]++;
	}

	/** increment age */
	for (unsigned k = 0; k < cols_generated_.size(); ++k) {
		if (cols_generated_[k]->active_) {
			int j = cols_generated_[k]->master_index_;
			if (getSiPtr()->getRowPrice()[j] < 1.0e-6)
				cols_generated_[k]->age_++;
			else
				cols_generated_[k]->age_ = 0;
			/** reduced cost fixing */
			if (cols_generated_[k]->age_ >= par_->getIntParam("DW/MASTER/COL_AGE_LIM") &&
				delblk[cols_generated_[k]->blockid_] > 1 &&
				getSiPtr()->getNumRows() - delrows.size() > getSiPtr()->getNumCols()) {

				cols_generated_[k]->active_ = false;
				delrows.push_back(j);
				delblk[cols_generated_[k]->blockid_]--;
			}
		}
	}

	if (delrows.size() > 0) {
		getSiPtr()->deleteRows(delrows.size(), &delrows[0]);
		primal_si_->si_->deleteCols(delrows.size(), &delrows[0]);

		/** reset age? */
		for (unsigned k = 0; k < cols_generated_.size(); ++k) {
			if (cols_generated_[k]->active_)
				cols_generated_[k]->age_ = 0;
		}
	}

	num_removed = delrows.size();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwBundleDual::restoreCols(int &num_restored) {

	num_restored = 0;
#if 0
	BGN_TRY_CATCH

	/** mark the existing columns to match with basis; and mark columns to delete */
	for (unsigned k = 0; k < cols_generated_.size(); ++k) {
		if (cols_generated_[k]->active_) continue;
		if (cols_generated_[k]->age_ < COIN_INT_MAX) {
			/** update column */
			DSP_RTN_CHECK_RTN_CODE(updateCol(cols_generated_[k]));

			/** add column */
			double lhs = cols_generated_[k]->col_.dotProduct(&dualsol_[0]);
			if (lhs - cols_generated_[k]->obj_ > 1.0e-6) {
				cols_generated_[k]->active_ = true;
				cols_generated_[k]->age_ = 0;
				/** set master index */
				cols_generated_[k]->master_index_ = primal_si_->si_->getNumCols();
				/** add row and column */
				addDualRow(cols_generated_[k]->col_, -COIN_DBL_MAX, cols_generated_[k]->obj_);
				primal_si_->si_->addCol(cols_generated_[k]->col_, 0.0, COIN_DBL_MAX, cols_generated_[k]->obj_);
				num_restored++;
			}
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
#endif
	return DSP_RTN_OK;
}

void DwBundleDual::removeAllCols() {
	/** remove all columns in the primal master */
	removeAllPrimCols();
	/** remove all rows in the dual master */
	removeAllDualRows();
}

void DwBundleDual::removeAllPrimCols() {
	std::vector<int> delcols(primal_si_->si_->getNumCols());
	std::iota(delcols.begin(), delcols.end(), 0.0);
	primal_si_->si_->deleteCols(primal_si_->si_->getNumCols(), &delcols[0]);
}

void DwBundleDual::removeAllDualRows() {
	std::vector<int> delrows(getSiPtr()->getNumRows());
	std::iota(delrows.begin(), delrows.end(), 0.0);
	getSiPtr()->deleteRows(getSiPtr()->getNumRows(), &delrows[0]);
}

void DwBundleDual::removeBranchingRows() {
	if (nrows_branch_ > 0) {
		std::vector<int> delinds(nrows_branch_);
		std::iota(delinds.begin(), delinds.end(), nrows_core_);
		primal_si_->si_->deleteRows(nrows_branch_, &delinds[0]);
		getSiPtr()->deleteCols(nrows_branch_, &delinds[0]);
		nrows_branch_ = 0;
	}
}

void DwBundleDual::addBranchingRow(double lb, double ub) {
	primal_si_->si_->addRow(0, NULL, NULL, lb, ub);
	getSiPtr()->addCol(0, NULL, NULL, 0.0, ub, lb);
}

void DwBundleDual::addBranchingCol(const CoinPackedVector& col, double obj) {
	addDualRow(col, -COIN_DBL_MAX, obj);
	primal_si_->si_->addCol(col, 0.0, COIN_DBL_MAX, obj);
}

void DwBundleDual::printIterInfo() {
	message_->print(1, "Iteration %3d: DW Bound %+e, ", itercnt_, primobj_);
	message_->print(2, "Dual %+e, Approx %+e, ", -dualobj_, -bestdualobj_-v_);
	message_->print(1, "Best Dual %+e ", -bestdualobj_);
	if (getRelApproxGap() < 1000)
		message_->print(1, "(gap %.2f %%), ", getRelApproxGap()*100);
	else
		message_->print(1, "(gap Large %%), ");
	message_->print(1, "nrows %d, ncols %d, ", getSiPtr()->getNumRows(), getSiPtr()->getNumCols());
	message_->print(1, "timing (total %.2f, master %.2f, gencols %.2f), statue %d\n",
			t_total_, t_master_, t_colgen_, status_);
	message_->print(2, "  predicted ascent %+e, |p| %+e, alpha %+e, linerr %+e, eps %+e, u %+e, counter %d\n",
			-v_, absp_, alpha_, -linerr_, eps_, u_, counter_);

	/** log */
	log_time_.push_back(CoinWallclockTime());
	log_bestdual_bounds_.push_back(-bestdualobj_);
	log_bestprim_bounds_.push_back(bestprimobj_);
}
