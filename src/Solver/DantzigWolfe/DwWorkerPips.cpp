/*
 * DwWorkerPips.cpp
 *
 *  Created on: Dec 15, 2016
 *      Author: Kibaek Kim
 */

// #define DSP_DEBUG
#include "Model/TssModel.h"
#include "Solver/DantzigWolfe/DwWorkerPips.h"
#include "Utility/DspMpi.h"
#include "PIPSIpmInterface.h"
#include "sFactoryAug.h"
// #include "sFactoryAugSchurLeaf.h"
#include "MehrotraStochSolver.h"

//#define SCALE_TAU

/** display option for PIPS */
int gOoqpPrintLevel = 0;

DwWorkerPips::DwWorkerPips(
		DecModel * model,
		DspParams * par,
		DspMessage * message,
		MPI_Comm comm):
DwWorkerMpi(model, par, message, comm), pips_(NULL) {
	tss_ = dynamic_cast<TssModel*>(model_);
}

DSP_RTN_CODE DwWorkerPips::receiver() {
	int signal;
	std::vector<int> indices;
	std::vector<int> statuses;
	std::vector<double> cxs;
	std::vector<double> objs;
	std::vector<CoinPackedVector*> sols;
	bool terminate = false;
	CoinError::printErrors_ = true;
	std::vector<double> dummy;
	std::vector<DwCol*> cols;

	BGN_TRY_CATCH

	while (terminate == false) {
		/** receive a signal */
		MPI_Bcast(&signal, 1, MPI_INT, 0, comm_);
		DSPdebugMessage("Rank %d received signal %d.\n", comm_rank_, signal);

		switch(signal) {
		case sig_generateCols:
			DSP_RTN_CHECK_RTN_CODE(
					generateCols(-1, NULL, indices, statuses, cxs, objs, sols));
			indices.clear();
			statuses.clear();
			cxs.clear();
			objs.clear();
			for (size_t i = 0; i < sols.size(); ++i) {
				delete sols[i];
				sols[i] = NULL;
			}
			sols.clear();
			break;
		case sig_sync:
			sync(dummy, cols);
			break;
		case sig_solvePips:
			DSP_RTN_CHECK_RTN_CODE(solvePips(0.0));
			break;
		case sig_setColBounds:
			setColBounds(0, NULL, NULL, NULL);
			break;
		case sig_terminate:
			terminate = true;
			break;
		default:
			CoinError("Unexpected signal", "DwWorkerPips", "coordinator");
			break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwWorkerPips::sync(std::vector<double> bestdualsol, std::vector<DwCol*> dwcols) {

	int num_blocks;
	int total_length;
	int num_elems;
	std::vector<int> block_ids;
	std::vector<int> x_lengths;
	std::vector<int> x_indices;
	std::vector<double> x_elements;
	std::vector<double> objs;

	BGN_TRY_CATCH

	if (tss_ == NULL) throw "tss_ is not defined.\n";

	if (comm_rank_ == 0) {
		/** send signal */
		int sig = sig_sync;
		DSPdebugMessage("Rank %d sends signal %d.\n", comm_rank_, sig);
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
	}

	if (comm_rank_ == 0) {
		num_blocks = 0;
		total_length = 0;
		block_ids.reserve(dwcols.size());
		x_lengths.reserve(dwcols.size());
		x_indices.reserve(dwcols.size() * tss_->getNumCols(0));
		x_elements.reserve(dwcols.size() * tss_->getNumCols(0));
		objs.reserve(dwcols.size());

		/** store column data */
		for (auto col = dwcols.begin(); col != dwcols.end(); col++)
			if ((*col)->active_) {
				num_blocks++;
				num_elems = 0;
				block_ids.push_back((*col)->blockid_);
				for (int j = 0; j < (*col)->x_.getNumElements(); ++j)
					if ((*col)->x_.getIndices()[j] < tss_->getNumCols(0) * tss_->getNumScenarios()) {
						x_indices.push_back((*col)->x_.getIndices()[j] % tss_->getNumCols(0));
						x_elements.push_back((*col)->x_.getElements()[j]);
						num_elems++;
					}
				x_lengths.push_back(num_elems);
				total_length += num_elems;
				objs.push_back((*col)->obj_);
			}
	}

	/** sync data */
	MPI_Bcast(&num_blocks, 1, MPI_INT, 0, comm_);
	MPI_Bcast(&total_length, 1, MPI_INT, 0, comm_);
	if (comm_rank_ > 0) {
		block_ids.resize(num_blocks);
		x_lengths.resize(num_blocks);
		objs.reserve(num_blocks);
		x_indices.resize(total_length);
		x_elements.resize(total_length);
		bestdualsol.resize((tss_->getNumCols(0) + 1) * tss_->getNumScenarios(), 0.0);
	}
	MPI_Bcast(&block_ids[0], num_blocks, MPI_INT, 0, comm_);
	MPI_Bcast(&x_lengths[0], num_blocks, MPI_INT, 0, comm_);
	MPI_Bcast(&x_indices[0], total_length, MPI_INT, 0, comm_);
	MPI_Bcast(&x_elements[0], total_length, MPI_DOUBLE, 0, comm_);
	MPI_Bcast(&objs[0], num_blocks, MPI_DOUBLE, 0, comm_);
	MPI_Bcast(&bestdualsol[0], (tss_->getNumCols(0) + 1) * tss_->getNumScenarios(), MPI_DOUBLE, 0, comm_);

	/** free the previous pips_ */
	FREE_PTR(pips_);

	/** create PIPS input data */
	pips_ = new DwBundlePipsInput(tss_->getNumScenarios(), tss_->getNumCols(0));

#if 0
	for (int i = 0; i < comm_size_; ++i) {
		if (i == comm_rank_) {
			printf("Rank %d\n", i);
			printf("Block ID:\n"); DspMessage::printArray(block_ids.size(), &block_ids[0]);
			printf("Length of x:\n"); DspMessage::printArray(x_lengths.size(), &x_lengths[0]);
			printf("Indices:\n"); DspMessage::printArray(x_indices.size(), &x_indices[0]);
			printf("Elements:\n"); DspMessage::printArray(x_elements.size(), &x_elements[0]);
		}
		MPI_Barrier(comm_);
	}

	if (comm_rank_ == 0) {
		printf("bestdualsol:\n");
		DspMessage::printArray(bestdualsol.size(), &bestdualsol[0]);
	}
#endif

	/** add proximal center */
	std::vector<std::vector<double>> prox_center(tss_->getNumScenarios());
	for (int s = 0; s < tss_->getNumScenarios(); ++s) {
		prox_center[s].assign(tss_->getNumCols(0), 0.0);
		for (int j = 0; j < tss_->getNumCols(0); ++j) {
			prox_center[s][j] = bestdualsol[tss_->getNumScenarios() + s * tss_->getNumCols(0) + j];
		}
	}
	pips_->setProxCenter(prox_center);

	/** add bundle information */
	int pos = 0;
	for (size_t s = 0; s < block_ids.size(); ++s) {
		// if (block_ids[s] % comm_size_ == comm_rank_)
			pips_->addBundleInfo(block_ids[s], x_lengths[s], &x_indices[pos], &x_elements[pos], objs[s]);
		pos += x_lengths[s];
	}

	// if (comm_rank_ == 0)
	// 	pips_->displayBundle();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
	return DSP_RTN_OK;
}

DSP_RTN_CODE DwWorkerPips::solvePips(double weight) {

	BGN_TRY_CATCH

	if (comm_rank_ == 0) {
		/** send signal */
		int sig = sig_solvePips;
		DSPdebugMessage("Rank %d sends signal %d.\n", comm_rank_, sig);
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
	}

	/** set tau */
	MPI_Bcast(&weight, 1, MPI_DOUBLE, 0, comm_);
	pips_->setTau(weight);

	/** Create PIPS-IPM object */
	if (comm_rank_ == 0) printf("creating PIPS-IPM object...\n");
	PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> solver(*pips_);
	// PIPSIpmInterface<sFactoryAugSchurLeaf, MehrotraStochSolver> solver(*input_);
fflush(stdout);
MPI_Barrier(comm_);

	if (comm_rank_ == 0) printf("solving the master with PIPS-IPM... \n");
	solver.go();
	pips_objval_ = solver.getObjective();
#ifdef SCALE_TAU
	pips_objval_ /= weight;
#endif
	if (comm_rank_ == 0) printf("done with objective %e\n", pips_objval_);

	/** store first-stage variable solution */
	std::vector<double> const& soln_w = solver.getFirstStagePrimalColSolution();
	assert(soln_w.size() == pips_->nvars1_);
#ifdef SCALE_TAU
	for (size_t j = 0; j < soln_w.size(); ++j)
		soln_w[j] /= weight;
#endif
	// if (comm_rank_ == 0) {
	// 	printf("w:\n"); DspMessage::printArray(soln_w.size(), &soln_w[0]);
	// }

	/** store second-stage variable solution */
	std::vector<std::vector<double>> soln_z(pips_->nscen_);
	std::vector<std::vector<double>> soln_theta(pips_->nscen_);
	for (int s = 0; s < pips_->nscen_; ++s) {
		soln_z[s] = solver.getSecondStagePrimalColSolution(s);
		soln_theta[s] = solver.getSecondStageDualRowSolution(s);
#ifdef SCALE_TAU
		for (size_t j = 0; j < soln_theta[s].size(); ++j)
			soln_theta[s][j] /= weight;
#endif
		// if (soln_z[s].size() > 0) {
		// 	printf("z[%d]:\n", s); 
		// 	DspMessage::printArray(soln_z[s].size(), &soln_z[s][0]);
		// }
	}

	/** calculate beta for each scenario */
	std::vector<std::vector<double>> soln_beta(pips_->nscen_);
	for (int s = 0; s < pips_->nscen_; ++s) {
		if (soln_z[s].size() > 0) {
			soln_beta[s].assign(tss_->getNumCols(0), 0.0);
			CoinPackedMatrix H = pips_->getSecondStageCrossHessian(s);
			assert(H.getNumRows() == soln_z[s].size());
			assert(H.getNumCols() == soln_beta[s].size());
			H.transposeTimes(&soln_z[s][0], &soln_beta[s][0]);
			assert(soln_w.size() <= soln_beta[s].size());
			for (size_t j = 0; j < soln_beta[s].size(); ++j) {
				// soln_beta[s][j] *= -1.0;
				soln_beta[s][j] += soln_w[j] / weight;
			}
		}
	}

	/** collect the second-stage solutions to root */
	int length_of_beta, length_of_theta;
	for (int s = 0; s < pips_->nscen_; ++s) {
		if (comm_rank_ == 0 && soln_beta[s].size() == 0) {
			MPI_Recv(&length_of_beta, 1, MPI_INT, MPI_ANY_SOURCE, 0, comm_, MPI_STATUS_IGNORE);
			MPI_Recv(&length_of_theta, 1, MPI_INT, MPI_ANY_SOURCE, 0, comm_, MPI_STATUS_IGNORE);
			soln_beta[s].resize(length_of_beta, 0.0);
			soln_theta[s].resize(length_of_theta, 0.0);
			MPI_Recv(&soln_beta[s][0], length_of_beta, MPI_DOUBLE, MPI_ANY_SOURCE, 0, comm_, MPI_STATUS_IGNORE);
			MPI_Recv(&soln_theta[s][0], length_of_theta, MPI_DOUBLE, MPI_ANY_SOURCE, 0, comm_, MPI_STATUS_IGNORE);
		} else if (comm_rank_ > 0 && soln_beta[s].size() > 0) {
			length_of_beta = soln_beta[s].size();
			length_of_theta = soln_theta[s].size();
			MPI_Send(&length_of_beta, 1, MPI_INT, 0, 0, comm_);
			MPI_Send(&length_of_theta, 1, MPI_INT, 0, 0, comm_);
			MPI_Send(&soln_beta[s][0], length_of_beta, MPI_DOUBLE, 0, 0, comm_);
			MPI_Send(&soln_theta[s][0], length_of_theta, MPI_DOUBLE, 0, 0, comm_);
		}
		MPI_Barrier(comm_);
	}

	if (comm_rank_ == 0) {
		pips_beta_ = soln_beta;
		pips_theta_ = soln_theta;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}
