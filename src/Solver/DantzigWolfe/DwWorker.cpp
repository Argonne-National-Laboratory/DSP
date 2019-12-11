/*
 * DwWorker.cpp
 *
 *  Created on: Aug 27, 2016
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

#ifdef DSP_HAS_CPX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#endif

#ifdef HAS_CBC
#include "OsiCbcSolverInterface.hpp"
#endif

#include "Solver/DantzigWolfe/DwWorker.h"
#include "Model/TssModel.h"
#include "Utility/DspUtility.h"

#define SI OsiCpxSolverInterface

DwWorker::DwWorker(DecModel * model, DspParams * par, DspMessage * message) :
		model_(model),
		par_(par),
		message_(message),
		si_(NULL),
		sub_objs_(NULL) {

	/** parameters */
	parProcIdxSize_ = par_->getIntPtrParamSize("ARR_PROC_IDX");
	parProcIdx_     = par_->getIntPtrParam("ARR_PROC_IDX");
	DSPdebugMessage("Created parameters, DwWorker.\n");

	if (parProcIdxSize_ <= 0) {
		throw printf("Parameter ARR_PROC_IDX should be assigned. in DwWorker::DwWorker\n");
	}

	/** number of total subproblems */
	nsubprobs_ = parProcIdxSize_;
	DSPdebugMessage("nsubprobs_ %d\n", nsubprobs_);

	/** create subproblem solver */
	//sub_ = new DwSub();

	/** create solver interface */
	si_ = new OsiSolverInterface* [parProcIdxSize_];
	for (int i = 0; i < parProcIdxSize_; ++i)
		si_[i] = new SI();

	/** subproblem objective coefficients */
	sub_objs_ = new double* [parProcIdxSize_];
	sub_clbd_ = new double* [parProcIdxSize_];
	sub_cubd_ = new double* [parProcIdxSize_];
	coupled_ = new bool* [parProcIdxSize_];
	for (int i = 0; i < parProcIdxSize_; ++i) {
		sub_objs_[i] = NULL;
		sub_clbd_[i] = NULL;
		sub_cubd_[i] = NULL;

		/** indicate whether columns are coupled with the master or not. */
		coupled_[i] = new bool [model_->getNumCouplingCols()];
		CoinFillN(coupled_[i], model_->getNumCouplingCols(), false);

		int nccols = model_->getNumSubproblemCouplingCols(parProcIdx_[i]);
		const int* ccols = model_->getSubproblemCouplingColIndices(parProcIdx_[i]);
		DSPdebugMessage("Subproblem(%d) coupling columns:\n", parProcIdx_[i]);
		DSPdebug(DspMessage::printArray(nccols, ccols));
		for (int j = 0; j < nccols; ++j)
			coupled_[i][ccols[j]] = true;
	}

	/** create subproblems */
	DSP_RTN_CHECK_THROW(createSubproblems());

	added_rowids_.resize(parProcIdxSize_);
}

DwWorker::~DwWorker() {
	FREE_2D_PTR(parProcIdxSize_, si_);
	//FREE_PTR(sub_);
	FREE_2D_ARRAY_PTR(parProcIdxSize_, sub_objs_);
	FREE_2D_ARRAY_PTR(parProcIdxSize_, sub_clbd_);
	FREE_2D_ARRAY_PTR(parProcIdxSize_, sub_cubd_);
	FREE_2D_ARRAY_PTR(parProcIdxSize_, coupled_);
}

DSP_RTN_CODE DwWorker::createSubproblems() {
#define FREE_MEMORY        \
	FREE_PTR(mat);    \
	FREE_ARRAY_PTR(ctype); \
	FREE_ARRAY_PTR(rlbd);  \
	FREE_ARRAY_PTR(rubd);

	TssModel* tss = NULL;
	CoinPackedMatrix* mat = NULL;
	char* ctype = NULL;
	double* rlbd = NULL;
	double* rubd = NULL;

	BGN_TRY_CATCH

	if (model_->isStochastic())
		tss = dynamic_cast<TssModel*>(model_);

	num_timelim_stops_.assign(parProcIdxSize_, 0);

	for (int s = 0; s < parProcIdxSize_; ++s) {
		if (model_->isStochastic()) {
			DSP_RTN_CHECK_RTN_CODE(
					model_->decompose(1, &parProcIdx_[s], 0, NULL, NULL, NULL,
							mat, sub_clbd_[s], sub_cubd_[s], ctype, sub_objs_[s], rlbd, rubd));
			for (int j = 0; j < tss->getNumCols(0); ++j)
				sub_objs_[s][j] *= tss->getProbability()[parProcIdx_[s]];
		} else {
			DSP_RTN_CHECK_RTN_CODE(
					model_->copySubprob(parProcIdx_[s], mat, sub_clbd_[s], sub_cubd_[s], ctype, sub_objs_[s], rlbd, rubd));
			DSPdebug(mat->verifyMtx(4));

			/** fix zeros for non-coupling columns */
			for (int j = 0; j < model_->getNumCouplingCols(); ++j) {
				if (coupled_[s][j] == false) {
					sub_clbd_[s][j] = 0.0;
					sub_cubd_[s][j] = 0.0;
					sub_objs_[s][j] = 0.0;
				}
			}
		}
		DSPdebugMessage("sub_objs_[%d]:\n", parProcIdx_[s]);
		DSPdebug(DspMessage::printArray(model_->getNumCouplingCols(), sub_objs_[s]));

		/** load problem to si */
		si_[s]->loadProblem(*mat, sub_clbd_[s], sub_cubd_[s], sub_objs_[s], rlbd, rubd);

		/** set integers */
		int nintegers = 0;
		for (int j = 0; j < si_[s]->getNumCols(); ++j) {
			if (ctype[j] != 'C') {
				si_[s]->setInteger(j);
				nintegers++;
			}
		}

		si_[s]->messageHandler()->setLogLevel(par_->getIntParam("DW/SUB/LOG_LEVEL"));

		/** set parameters */
		setGapTolerance(par_->getDblParam("DW/SUB/GAPTOL"));
		setTimeLimit(par_->getDblParam("DW/SUB/TIME_LIM"));
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_[s]);
		if (cpx) {
			if (nintegers > 0)
				cpx->switchToMIP();
			CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_THREADS, par_->getIntParam("DW/SUB/THREADS"));
			CPXsetintparam(cpx->getEnvironmentPtr(), CPX_PARAM_ADVIND, par_->getIntParam("DW/SUB/ADVIND"));
		}

		if (nintegers > 0)
			si_[s]->branchAndBound();
		else
			si_[s]->initialSolve();

#ifdef DSP_DEBUG
		if (s >= 0) {
			/** write MPS */
			char ofname[128];
			sprintf(ofname, "sub%d", parProcIdx_[s]);
			DSPdebugMessage("Writing MPS file: %s\n", ofname);
			si_[s]->writeMps(ofname);
		}
#endif

		/** free memory */
		FREE_MEMORY
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** generate variables */
DSP_RTN_CODE DwWorker::generateCols(
		int phase,                           /**< [in] phase of the master */
		const double* piA,                   /**< [in] piA */
		std::vector<int>& indices,           /**< [out] subproblem indices */
		std::vector<int>& statuses,          /**< [out] solution status */
		std::vector<double>& cxs,            /**< [out] solution times original objective coefficients */
		std::vector<double>& objs,           /**< [out] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */) {
	CoinError::printErrors_ = true;

	/** subproblem objective and solution */
	double cx;
	double objval;
	CoinPackedVector* sol = NULL;
	TssModel* tss = NULL;

	BGN_TRY_CATCH

	/** adjust objective function */
	DSP_RTN_CHECK_RTN_CODE(adjustObjFunction(phase, piA));

	/** solve subproblems */
	DSP_RTN_CHECK_RTN_CODE(solveSubproblems());

	/** cleanup and reserve memory*/
	indices.clear();
	statuses.clear();
	cxs.clear();
	objs.clear();
	for (unsigned i = 0; i < sols.size(); ++i)
		FREE_PTR(sols[i]);
	sols.clear();
	indices.reserve(parProcIdxSize_);
	statuses.reserve(parProcIdxSize_);
	cxs.reserve(parProcIdxSize_);
	objs.reserve(parProcIdxSize_);
	sols.reserve(parProcIdxSize_);


	if (model_->isStochastic())
		tss = dynamic_cast<TssModel*>(model_);

	for (int s = 0; s < parProcIdxSize_; ++s) {
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_[s]);
		int sind = parProcIdx_[s];

		/** add subproblem index */
		indices.push_back(sind);

		/** store solution status */
		int status;
		int cpxstat = CPXgetstat(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
		convertCoinToDspStatus(si_[s], status);
		/** FIXME: Osi does not know this. */
		if (cpxstat == CPXMIP_INFEASIBLE) {
			message_->print(2, "  Subproblem %d infeasible.\n", sind);
			status = DSP_STAT_PRIM_INFEASIBLE;
		} else if (cpxstat == CPXMIP_INForUNBD) {
			message_->print(2, "  Subproblem %d unbounded.\n", sind);
			status = DSP_STAT_DUAL_INFEASIBLE;
		} else if (cpxstat == CPXMIP_TIME_LIM_FEAS) {
			message_->print(2, "  Subproblem %d terminated due to time limit.\n", sind);
			status = DSP_STAT_LIM_ITERorTIME;
			num_timelim_stops_[s]++;
		} else if (status != DSP_STAT_UNKNOWN) {
			num_timelim_stops_[s] = 0;
			setTimeLimit(par_->getDblParam("DW/SUB/TIME_LIM"));
		}

		DSPdebugMessage("sind %d status %d\n", sind, status);
		statuses.push_back(status);

		if (status != DSP_STAT_UNKNOWN && status != DSP_STAT_PRIM_INFEASIBLE) {
			sol = new CoinPackedVector;
			sol->reserve(si_[s]->getNumCols());

			if (status == DSP_STAT_DUAL_INFEASIBLE) {
				/** retrieve ray if unbounded */
				std::vector<double*> rays = si_[s]->getPrimalRays(1);
				if (rays.size() == 0 || rays[0] == NULL)
					throw CoinError("No primal ray is available.", "generateCols", "DwWorker");
				double* ray = rays[0];
				rays[0] = NULL;

				/** subproblem objective value */
				cx = 0.0;
				objval = 0.0;
				for (int j = 0; j < si_[s]->getNumCols(); ++j) {
					cx += sub_objs_[s][j] * ray[j];
					objval += si_[s]->getObjCoefficients()[j] * ray[j];
				}

				/** subproblem coupling solution */
				for (int j = 0; j < si_[s]->getNumCols(); ++j) {
					double xval = ray[j];
					if (fabs(xval) > 1.0e-8) {
						if (model_->isStochastic()) {
							if (j < tss->getNumCols(0))
								sol->insert(sind * tss->getNumCols(0) + j, xval);
							else
								sol->insert((tss->getNumScenarios()-1) * tss->getNumCols(0) + sind * tss->getNumCols(1) + j, xval);
						} else
							sol->insert(j, xval);
					}
				}

				/** free ray */
				FREE_ARRAY_PTR(ray);
			} else /*if (si_[s]->isProvenOptimal())*/ {
				const double* x = si_[s]->getColSolution();

				/** subproblem objective value */
				if (si_[s]->getNumIntegers() > 0)
					CPXgetbestobjval(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &objval);
				else
					objval = si_[s]->getObjValue();
					
				cx = 0.0;
				for (int j = 0; j < si_[s]->getNumCols(); ++j)
					cx += sub_objs_[s][j] * x[j];
				DSPdebugMessage("Subprob %d: objval %e, cx %e\n", sind, objval, cx);

				/** subproblem coupling solution */
				for (int j = 0; j < si_[s]->getNumCols(); ++j) {
					double xval = x[j];
					//if (sub_clbd_[s][j] == sub_cubd_[s][j])
					//	printf("sind %d j %d [%e, %e, %e]\n", sind, j, sub_clbd_[s][j], xval, sub_cubd_[s][j]);
					if (fabs(xval) > 1.0e-8) {
						if (model_->isStochastic()) {
							if (j < tss->getNumCols(0))
								sol->insert(sind * tss->getNumCols(0) + j, xval);
							else
								sol->insert((tss->getNumScenarios()-1) * tss->getNumCols(0) + sind * tss->getNumCols(1) + j, xval);
						} else
							sol->insert(j, xval);
					}
				}
			}

			/** store objective and solution */
			cxs.push_back(cx);
			objs.push_back(objval);
			sols.push_back(sol);
			sol = NULL;
		} else {
			message_->print(1, "Unexpected subproblem status (block: %d, status: %d)\n", sind, status);
			/** store dummies */
			cxs.push_back(0.0);
			objs.push_back(0.0);
			sols.push_back(new CoinPackedVector);
		}
	}

	/** reset subproblems */
	DSP_RTN_CHECK_RTN_CODE(resetSubproblems());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwWorker::generateColsByFix(
		const double* x,                     /**< [in] solution to fix */
		std::vector<int>& indices,           /**< [out] subproblem indices */
		std::vector<int>& statuses,          /**< [out] solution status */
		std::vector<double>& objs,           /**< [out] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */) {

	/** run only for stochastic models */
	if (model_->isStochastic() == false)
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	double objval;
	CoinPackedVector* sol = NULL;
	TssModel* tss = dynamic_cast<TssModel*>(model_);
	int ncols_first_stage = tss->getNumCols(0);

	/** set objective function */
	for (int s = 0; s < parProcIdxSize_; ++s) {
		for (int j = 0; j < si_[s]->getNumCols(); ++j) {
			si_[s]->setObjCoeff(j, sub_objs_[s][j]);
			// message_->print(0, "sub_objs_[%d][%d] = %e\n", s, j, sub_objs_[s][j]);
		}
	}

	/** fix column bounds */
	for (int s = 0; s < parProcIdxSize_; ++s) {
		for (int j = 0; j < ncols_first_stage; ++j)
			si_[s]->setColBounds(j, x[j], x[j]);
	}

	/** solve subproblems */
	DSP_RTN_CHECK_RTN_CODE(solveSubproblems());

	/** cleanup and reserve memory*/
	indices.clear();
	statuses.clear();
	objs.clear();
	for (unsigned i = 0; i < sols.size(); ++i)
		FREE_PTR(sols[i]);
	sols.clear();
	indices.reserve(parProcIdxSize_);
	statuses.reserve(parProcIdxSize_);
	objs.reserve(parProcIdxSize_);
	sols.reserve(parProcIdxSize_);

	for (int s = 0; s < parProcIdxSize_; ++s) {
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_[s]);
		int sind = parProcIdx_[s];

		/** add subproblem index */
		indices.push_back(sind);

		/** store solution status */
		int status;
		int cpxstat = CPXgetstat(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
		convertCoinToDspStatus(si_[s], status);
		/** FIXME: Osi does not know this. */
		if (cpxstat == CPXMIP_INFEASIBLE) {
			message_->print(3, "  Upper bounding subproblem %d infeasible.\n", sind);
			status = DSP_STAT_PRIM_INFEASIBLE;
		} else if (cpxstat == CPXMIP_INForUNBD) {
			message_->print(3, "  Upper bounding subproblem %d unbounded.\n", sind);
			status = DSP_STAT_DUAL_INFEASIBLE;
		} else if (cpxstat == CPXMIP_TIME_LIM_FEAS) {
			message_->print(3, "  Upper bounding subproblem %d terminated due to time limit.\n", sind);
			status = DSP_STAT_LIM_ITERorTIME;
			num_timelim_stops_[s]++;
		} else if (status != DSP_STAT_UNKNOWN) {
			num_timelim_stops_[s] = 0;
			setTimeLimit(par_->getDblParam("DW/SUB/TIME_LIM"));
		}

		DSPdebugMessage("sind %d status %d\n", sind, status);
		statuses.push_back(status);

		if (!si_[s]->isAbandoned() && !si_[s]->isProvenPrimalInfeasible()) {
			sol = new CoinPackedVector;
			sol->reserve(si_[s]->getNumCols());

			if (!si_[s]->isProvenDualInfeasible()) {
				const double* x = si_[s]->getColSolution();

				/** subproblem objective value */
				if (si_[s]->getNumIntegers())
					CPXgetbestobjval(cpx->getEnvironmentPtr(), cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &objval);
				else
					objval = si_[s]->getObjValue();
				DSPdebugMessage("Subprob %d: objval %e\n", sind, objval);

				/** subproblem coupling solution */
				for (int j = 0; j < si_[s]->getNumCols(); ++j) {
					double xval = x[j];
					//if (sub_clbd_[s][j] == sub_cubd_[s][j])
					//	printf("sind %d j %d [%e, %e, %e]\n", sind, j, sub_clbd_[s][j], xval, sub_cubd_[s][j]);
					if (fabs(xval) > 1.0e-8) {
						if (j < tss->getNumCols(0))
							sol->insert(sind * tss->getNumCols(0) + j, xval);
						else
							sol->insert((tss->getNumScenarios()-1) * tss->getNumCols(0) + sind * tss->getNumCols(1) + j, xval);
					}
				}
			}

			/** store objective and solution */
			objs.push_back(objval);
			sols.push_back(sol);
			sol = NULL;
		} else {
			message_->print(0, "Unexpected subproblem status (block: %d, status: %d)\n", sind, status);
			/** store dummies */
			objs.push_back(0.0);
			sols.push_back(new CoinPackedVector);
		}
	}

	/** reset subproblems */
	DSP_RTN_CHECK_RTN_CODE(resetSubproblems());

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwWorker::adjustObjFunction(
		int phase,        /**< [in] phase of the master */
		const double* piA /**< [in] dual variable times the constraint matrix */) {

	TssModel* tss = NULL;
	int ncols_first_stage, ncols_second_stage, nscen;

	BGN_TRY_CATCH

	if (model_->isStochastic()) {
		tss = dynamic_cast<TssModel*>(model_);
		ncols_first_stage = tss->getNumCols(0);
		ncols_second_stage = tss->getNumCols(1);
		nscen = tss->getNumScenarios();
	}

	DSPdebugMessage("adjustObjFunction is in phase %d.\n", phase);
	for (int s = 0; s < parProcIdxSize_; ++s) {
		/** actual subproblem index */
		int sind = parProcIdx_[s];

		/** set new objective coefficients */
		if (model_->isStochastic()) {
			if (phase == 1) {
				for (int j = 0; j < ncols_first_stage; ++j)
					si_[s]->setObjCoeff(j, -piA[sind * ncols_first_stage + j]);
				for (int j = ncols_first_stage; j < si_[s]->getNumCols(); ++j)
					si_[s]->setObjCoeff(j, -piA[(nscen-1) * ncols_first_stage + sind * ncols_second_stage + j]);
			} else if (phase == 2) {
				for (int j = 0; j < ncols_first_stage; ++j)
					si_[s]->setObjCoeff(j, sub_objs_[s][j] - piA[sind * ncols_first_stage + j]);
				for (int j = ncols_first_stage; j < si_[s]->getNumCols(); ++j)
					si_[s]->setObjCoeff(j, sub_objs_[s][j] - piA[(nscen-1) * ncols_first_stage + sind * ncols_second_stage + j]);
			}
		} else {
			if (phase == 1) {
				for (int j = 0; j < si_[s]->getNumCols(); ++j) {
					if (j < model_->getNumCouplingCols())
						si_[s]->setObjCoeff(j, -piA[j]);
					else
						si_[s]->setObjCoeff(j, 0.0);
				}
			} else if (phase == 2) {
				for (int j = 0; j < si_[s]->getNumCols(); ++j) {
					if (j < model_->getNumCouplingCols())
						si_[s]->setObjCoeff(j, sub_objs_[s][j] - piA[j]);
					else
						si_[s]->setObjCoeff(j, sub_objs_[s][j]);
				}
			}
		}
		//DSPdebugMessage("[Phase %d] Set objective coefficients for subproblem %d:\n", phase, sind);
		//DSPdebug(DspMessage::printArray(si_[s]->getNumCols(), si_[s]->getObjCoefficients()));
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** set column lower bounds */
void DwWorker::setColBounds(int j, double lb, double ub) {
	if (model_->isStochastic()) {
		TssModel* tss = dynamic_cast<TssModel*>(model_);
		for (int s = 0; s < parProcIdxSize_; ++s) {
			if (j < tss->getNumScenarios() * tss->getNumCols(0)) {
				si_[s]->setColBounds(j % tss->getNumCols(0), lb, ub);
				sub_clbd_[s][j % tss->getNumCols(0)] = lb;
				sub_cubd_[s][j % tss->getNumCols(0)] = ub;
				DSPdebugMessage("subproblem %d changed column bounds: %d [%e %e]\n", parProcIdx_[s], j % tss->getNumCols(0), lb, ub);
			} else {
				int jj = j - tss->getNumScenarios() * tss->getNumCols(0);
				if (jj < parProcIdx_[s] * tss->getNumCols(1) || jj >= (parProcIdx_[s]+1) * tss->getNumCols(1))
					continue;
				si_[s]->setColBounds(tss->getNumCols(0) + jj % tss->getNumCols(1), lb, ub);
				sub_clbd_[s][tss->getNumCols(0) + jj % tss->getNumCols(1)] = lb;
				sub_cubd_[s][tss->getNumCols(0) + jj % tss->getNumCols(1)] = ub;
				DSPdebugMessage("subproblem %d changed column bounds: %d [%e %e]\n", parProcIdx_[s], tss->getNumCols(0) + jj % tss->getNumCols(1), lb, ub);
			}
		}
	} else {
		for (int s = 0; s < parProcIdxSize_; ++s) {
			if (coupled_[s][j] == true) {
				si_[s]->setColBounds(j, lb, ub);
				sub_clbd_[s][j] = lb;
				sub_cubd_[s][j] = ub;
				DSPdebugMessage("subproblem %d changed column bounds: %d [%e %e]\n", parProcIdx_[s], j, lb, ub);
			}
		}
	}
}

void DwWorker::setColBounds(int size, const int* indices, const double* lbs, const double* ubs) {
	for (int i = 0; i < size; ++i)
		setColBounds(indices[i], lbs[i], ubs[i]);
}

void DwWorker::addRow(const CoinPackedVector* vec, double lb, double ub) {
	for (int s = 0; s < parProcIdxSize_; ++s) {
		si_[s]->addRow(*vec, lb, ub);
		added_rowids_[s].push_back(si_[s]->getNumRows()-1);
		//printf("addRow: s %d size %u row index %d\n", s, added_rowids_[s].size(), si_[s]->getNumRows()-1);
	}
}

void DwWorker::removeAddedRows() {
	for (int s = 0; s < parProcIdxSize_; ++s) {
		//printf("removeAddedRows (s %u)\n", added_rowids_[s].size());
		if (added_rowids_[s].size() > 0) {
			//DspMessage::printArray(added_rowids_[s].size(), &(added_rowids_[s])[0]);
			si_[s]->deleteRows(added_rowids_[s].size(), &(added_rowids_[s])[0]);
			added_rowids_[s].clear();
		}
	}
}

DSP_RTN_CODE DwWorker::solveSubproblems() {
	int status;

	BGN_TRY_CATCH

#ifdef DSP_DEBUG
	for (int s = 0; s < parProcIdxSize_; ++s) {
		if (s >= 0) {
			/** write MPS */
			char ofname[128];
			sprintf(ofname, "sub%d", s);
			DSPdebugMessage("Writing MPS file: %s\n", ofname);
			si_[s]->writeMps(ofname);
		}
	}
#endif

	double timlim = par_->getDblParam("DW/SUB/TIME_LIM");
	int max_stops = *std::max_element(num_timelim_stops_.begin(), num_timelim_stops_.end());
	max_stops = 0;
	if (max_stops > 0) {
		timlim *= std::pow(2,max_stops);
		message_->print(3, "  Increased the time limit to %f for subproblems\n", timlim);
	}

	/** TODO: That's it? Dual infeasible??? */
	for (int s = 0; s < parProcIdxSize_; ++s) {
#ifdef HAS_CBC
		/** reset problem status */
		const OsiCbcSolverInterface* cbc = dynamic_cast<OsiCbcSolverInterface*>(si_[s]);
		if (cbc)
			cbc->getModelPtr()->setProblemStatus(-1);
#endif

		if (si_[s]->getNumIntegers() > 0) {

			/** increase time limit */
			if (max_stops > 0)
				setTimeLimit(timlim);

			/** solve */
			si_[s]->branchAndBound();
#ifdef DSP_DEBUG
			convertCoinToDspStatus(si_[s], status);
			DSPdebugMessage("MILP subproblem %d status %d\n", parProcIdx_[s], status);
#endif
			if (si_[s]->isProvenDualInfeasible()) {
				/** If primal unbounded, ray may not be immediately available.
				 * But, it becomes available if it is solved one more time.
				 * This is probably because the resolve() above behaved as initialSolve(),
				 * in which case presolve determines unboundedness without solve.
				 */
				si_[s]->resolve();
			}
		} else {
			/** solve LP relaxation */
			si_[s]->resolve();
#ifdef DSP_DEBUG
			convertCoinToDspStatus(si_[s], status);
			DSPdebugMessage("LP relaxation subproblem %d status %d\n", parProcIdx_[s], status);
#endif
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

void DwWorker::setTimeLimit(double limit) {
	for (int s = 0; s < parProcIdxSize_; ++s) {
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_[s]);
		if (cpx)
			CPXsetdblparam(cpx->getEnvironmentPtr(), CPX_PARAM_TILIM, limit);
	}
}

void DwWorker::setGapTolerance(double gaptol) {
	for (int s = 0; s < parProcIdxSize_; ++s) {
		OsiCpxSolverInterface* cpx = dynamic_cast<OsiCpxSolverInterface*>(si_[s]);
		if (cpx)
			CPXsetdblparam(cpx->getEnvironmentPtr(), CPX_PARAM_EPGAP, gaptol);
	}
}

void DwWorker::resetTimeIncrement() {
	std::fill(num_timelim_stops_.begin(), num_timelim_stops_.end(), 0);
}

DSP_RTN_CODE DwWorker::resetSubproblems() {
	BGN_TRY_CATCH
	/** restore bounds */
	for (int s = 0; s < parProcIdxSize_; ++s) {
		for (int j = 0; j < si_[s]->getNumCols(); ++j)
			si_[s]->setColBounds(j, sub_clbd_[s][j], sub_cubd_[s][j]);
	}
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}
