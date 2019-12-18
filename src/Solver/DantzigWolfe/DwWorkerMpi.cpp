/*
 * DwWorkerMpi.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG
#include "Model/TssModel.h"
#include "Solver/DantzigWolfe/DwWorkerMpi.h"

DwWorkerMpi::DwWorkerMpi(
		DecModel * model,
		DspParams * par,
		DspMessage * message,
		MPI_Comm comm):
DwWorker(model, par, message),
comm_(comm), resetTimeIncrement_(0) {
	MPI_Comm_rank(comm_, &comm_rank_);
	MPI_Comm_size(comm_, &comm_size_);

	/** calculate the size of piA vector */
	if (model_->isStochastic()) {
		TssModel* tss = dynamic_cast<TssModel*>(model_);
		npiA_ = tss->getNumScenarios() * (tss->getNumCols(0) + tss->getNumCols(1));
	} else
		npiA_ = model_->getNumCouplingCols();

	/** create a local vector memory */
	piA_ = new double [npiA_];
	DSPdebugMessage("Rank %d created piA_ vector of size %d.\n", comm_rank_, npiA_);

	/** get the number of subproblems */
	MPI_Allreduce(&parProcIdxSize_, &nsubprobs_, 1, MPI_INT, MPI_SUM, comm_);
#ifdef DSP_DEBUG
	if (comm_rank_ == 0)
		DSPdebugMessage("Number of subproblems: %d\n", nsubprobs_);
#endif
}

DwWorkerMpi::~DwWorkerMpi() {
	FREE_ARRAY_PTR(piA_);
}

DSP_RTN_CODE DwWorkerMpi::receiver() {
	int signal;
	std::vector<int> indices;
	std::vector<int> statuses;
	std::vector<double> cxs;
	std::vector<double> objs;
	std::vector<CoinPackedVector*> sols;
	bool terminate = false;
	CoinError::printErrors_ = true;

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
		case sig_generateColsByFix:
			DSP_RTN_CHECK_RTN_CODE(
					generateColsByFix(NULL, indices, statuses, objs, sols));
			indices.clear();
			statuses.clear();
			cxs.clear();
			objs.clear();
			sols.clear();
			break;
		case sig_setColBounds:
			setColBounds(0, NULL, NULL, NULL);
			break;
		case sig_addRow:
			addRow(NULL, 0.0, 0.0);
			break;
		case sig_terminate:
			terminate = true;
			break;
		default:
			CoinError("Unexpected signal", "DwWorkerMpi", "coordinator");
			break;
		}
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DwWorkerMpi::generateCols(
		int phase,                           /**< [in] phase of the master */
		const double* piA,                   /**< [in] piA */
		std::vector<int>& indices,           /**< [out] subproblem indices */
		std::vector<int>& statuses,          /**< [out] solution status */
		std::vector<double>& cxs,            /**< [out] solution times original objective coefficients */
		std::vector<double>& objs,           /**< [out] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */) {
#define FREE_MEMORY \
		FREE_ARRAY_PTR(recvcounts) \
		FREE_ARRAY_PTR(displs) \
		FREE_ARRAY_PTR(_indices) \
		FREE_ARRAY_PTR(_statuses) \
		FREE_ARRAY_PTR(_cxs) \
		FREE_ARRAY_PTR(_objs) \
		for (unsigned i = 0; i < _sols.size(); ++i) { \
			FREE_PTR(_sols[i]) \
		}

	int* recvcounts = NULL; /**< number of subproblems for each rank */
	int* displs = NULL;     /**< displacement for each receive buffer */

	/** temporary arrays to collect data from ranks */
	int* _indices = NULL;
	int* _statuses = NULL;
	double* _cxs = NULL;
	double* _objs = NULL;
	std::vector<CoinPackedVector*> _sols;

	BGN_TRY_CATCH

	if (comm_rank_ == 0) {
		/** send signal */
		int sig = sig_generateCols;
		DSPdebugMessage("Rank %d sends signal %d.\n", comm_rank_, sig);
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);

		/** copy piA vector */
		CoinCopyN(piA, npiA_, piA_);

		recvcounts = new int [comm_size_];
		displs = new int [comm_size_];
		_indices = new int [nsubprobs_];
		_statuses = new int [nsubprobs_];
		_cxs = new double [nsubprobs_];
		_objs = new double [nsubprobs_];
	}

	/** synchronize log level parameter */
	MPI_Bcast(&(message_->logLevel_), 1, MPI_INT, 0, comm_);

	/** synchronize phase and piA_ */
	MPI_Bcast(&phase, 1, MPI_INT, 0, comm_);
	MPI_Bcast(piA_, npiA_, MPI_DOUBLE, 0, comm_);
	DSPdebugMessage("Rank %d received phase %d and piA_ vector of size %d.\n", comm_rank_, phase, npiA_);

	/** reset time increment? */
	MPI_Bcast(&resetTimeIncrement_, 1, MPI_INT, 0, comm_);
	if (resetTimeIncrement_) {
		DwWorker::resetTimeIncrement();
		resetTimeIncrement_ = 0;
	}

	/** get maximum number of stops due to time limit */
	int maxstops_s = num_timelim_stops_.size() > 0 ? *std::max_element(num_timelim_stops_.begin(), num_timelim_stops_.end()) : 0;
	int maxstops_r = 0;
	DSPdebugMessage("Rank %d: maxstops_s %d maxstops_r %d\n", comm_rank_, maxstops_s, maxstops_r);
	MPI_Allreduce(&maxstops_s, &maxstops_r, 1, MPI_INT, MPI_MAX, comm_);
	std::fill(num_timelim_stops_.begin(), num_timelim_stops_.end(), maxstops_r);

	DSP_RTN_CHECK_RTN_CODE(
			DwWorker::generateCols(phase, piA_, indices, statuses, cxs, objs, sols));
	DSPdebugMessage("Rank %d generated %u indices, %u statuses, %u cxs, %u objs, and %u sols.\n",
			comm_rank_, indices.size(), statuses.size(), cxs.size(), objs.size(), sols.size());

	/** The root rank gathers the number of subproblems for each process. */
	MPI_Gather(&parProcIdxSize_, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, comm_);

	/** calculate displacement of the receive buffer */
	if (comm_rank_ == 0) {
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + recvcounts[i-1];
	}

	/** synchronize information */
	MPI_Gatherv(&indices[0], parProcIdxSize_, MPI_INT,
			_indices, recvcounts, displs, MPI_INT, 0, comm_);
	MPI_Gatherv(&statuses[0], parProcIdxSize_, MPI_INT,
			_statuses, recvcounts, displs, MPI_INT, 0, comm_);
	MPI_Gatherv(&cxs[0], parProcIdxSize_, MPI_DOUBLE,
			_cxs, recvcounts, displs, MPI_DOUBLE, 0, comm_);
	MPI_Gatherv(&objs[0], parProcIdxSize_, MPI_DOUBLE,
			_objs, recvcounts, displs, MPI_DOUBLE, 0, comm_);
	MPIgatherCoinPackedVectors(comm_, sols, _sols);

	/** construct output arguments */
	if (comm_rank_ == 0) {
		indices.resize(nsubprobs_);
		statuses.resize(nsubprobs_);
		cxs.resize(nsubprobs_);
		objs.resize(nsubprobs_);
		sols.resize(nsubprobs_);
		for (int i = 0; i < nsubprobs_; ++i) {
			indices[_indices[i]] = _indices[i];
			statuses[_indices[i]] = _statuses[i];
			cxs[_indices[i]] = _cxs[i];
			objs[_indices[i]] = _objs[i];
			DSPdebugMessage("_objs[%d] = %e\n", i, _objs[i]);
			sols[_indices[i]] = _sols[i];
			_sols[i] = NULL;
			// printf("i %d: indices %d, cxs %e, objs %e\n", i, _indices[i], _cxs[i], _objs[i]);
		}
		_sols.clear();
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

DSP_RTN_CODE DwWorkerMpi::generateColsByFix(
		const double* x,                     /**< [in] solution to fix */
		std::vector<int>& indices,           /**< [out] subproblem indices */
		std::vector<int>& statuses,          /**< [out] solution status */
		std::vector<double>& objs,           /**< [out] subproblem objective values */
		std::vector<CoinPackedVector*>& sols /**< [out] subproblem coupling column solutions */) {
#define FREE_MEMORY \
		FREE_ARRAY_PTR(recvcounts) \
		FREE_ARRAY_PTR(displs) \
		FREE_ARRAY_PTR(_x) \
		FREE_ARRAY_PTR(_indices) \
		FREE_ARRAY_PTR(_statuses) \
		FREE_ARRAY_PTR(_objs) \
		for (unsigned i = 0; i < _sols.size(); ++i) { \
			FREE_PTR(_sols[i]); \
		}

	int* recvcounts = NULL; /**< number of subproblems for each rank */
	int* displs = NULL;     /**< displacement for each receive buffer */

	/** temporary arrays to collect data from ranks */
	double* _x = NULL;
	int* _indices = NULL;
	int* _statuses = NULL;
	double* _objs = NULL;
	std::vector<CoinPackedVector*> _sols;

	/** run only for stochastic models */
	if (model_->isStochastic() == false)
		return DSP_RTN_OK;

	BGN_TRY_CATCH

	TssModel* tss = dynamic_cast<TssModel*>(model_);

	/** memory for solution to evalute */
	_x = new double [tss->getNumCols(0)];

	if (comm_rank_ == 0) {
		/** send signal */
		int sig = sig_generateColsByFix;
		DSPdebugMessage("Rank %d sends signal %d.\n", comm_rank_, sig);
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);

		recvcounts = new int [comm_size_];
		displs = new int [comm_size_];
		_indices = new int [nsubprobs_];
		_statuses = new int [nsubprobs_];
		_objs = new double [nsubprobs_];

		CoinCopyN(x, tss->getNumCols(0), _x);
	}

	MPI_Bcast(_x, tss->getNumCols(0), MPI_DOUBLE, 0, comm_);

	/** actual function to generate columns */
	DSP_RTN_CHECK_RTN_CODE(
			DwWorker::generateColsByFix(_x, indices, statuses, objs, sols));
	DSPdebugMessage("Rank %d generated %u indices, %u statuses, %u objs, and %u sols.\n",
			comm_rank_, indices.size(), statuses.size(), objs.size(), sols.size());

	/** The root rank gathers the number of subproblems for each process. */
	MPI_Gather(&parProcIdxSize_, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, comm_);

	/** calculate displacement of the receive buffer */
	if (comm_rank_ == 0) {
		displs[0] = 0;
		for (int i = 1; i < comm_size_; ++i)
			displs[i] = displs[i-1] + recvcounts[i-1];
	}

	/** synchronize information */
	MPI_Gatherv(&indices[0], parProcIdxSize_, MPI_INT,
			_indices, recvcounts, displs, MPI_INT, 0, comm_);
	MPI_Gatherv(&statuses[0], parProcIdxSize_, MPI_INT,
			_statuses, recvcounts, displs, MPI_INT, 0, comm_);
	MPI_Gatherv(&objs[0], parProcIdxSize_, MPI_DOUBLE,
			_objs, recvcounts, displs, MPI_DOUBLE, 0, comm_);
	MPIgatherCoinPackedVectors(comm_, sols, _sols);

	/** construct output arguments */
	if (comm_rank_ == 0) {
		indices.clear();
		statuses.clear();
		objs.clear();
		sols.clear();
		for (int i = 0; i < nsubprobs_; ++i) {
			indices.push_back(_indices[i]);
			statuses.push_back(_statuses[i]);
			objs.push_back(_objs[i]);
			sols.push_back(_sols[i]);
			_sols[i] = NULL;
		}
		_sols.clear();
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

void DwWorkerMpi::setColBounds(int size, const int* indices, const double* lbs, const double* ubs) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(_indices) \
	FREE_ARRAY_PTR(_lbs) \
	FREE_ARRAY_PTR(_ubs)

	int _size;
	int* _indices = NULL;
	double* _lbs = NULL;
	double* _ubs = NULL;

	BGN_TRY_CATCH

	/** send signal */
	if (comm_rank_ == 0) {
		int sig = sig_setColBounds;
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
		DSPdebugMessage("Rank 0 sent signal %d.\n", sig);
		_size = size;
		_indices = const_cast<int*>(indices);
		_lbs = const_cast<double*>(lbs);
		_ubs = const_cast<double*>(ubs);
	}

	/** synchronize information */
	MPI_Bcast(&_size, 1, MPI_INT, 0, comm_);
	DSPdebugMessage("Rank %d broadcasted size [%d].\n", comm_rank_, _size);
	if (comm_rank_ > 0) {
		_indices = new int [_size];
		_lbs = new double [_size];
		_ubs = new double [_size];
	}
	MPI_Bcast(_indices, _size, MPI_INT, 0, comm_);
	MPI_Bcast(_lbs, _size, MPI_DOUBLE, 0, comm_);
	MPI_Bcast(_ubs, _size, MPI_DOUBLE, 0, comm_);
	DwWorker::setColBounds(_size, _indices, _lbs, _ubs);

	if (comm_rank_ == 0) {
		_indices = NULL;
		_lbs = NULL;
		_ubs = NULL;
	}

	END_TRY_CATCH(FREE_MEMORY)
	FREE_MEMORY
#undef FREE_MEMORY
}

void DwWorkerMpi::addRow(const CoinPackedVector* vec, double lb, double ub) {
#define FREE_MEMORY \
	FREE_ARRAY_PTR(_indices) \
	FREE_ARRAY_PTR(_elements) \
	FREE_PTR(_vec)

	int _size;
	int* _indices = NULL;
	double* _elements = NULL;
	double _lb, _ub;
	CoinPackedVector* _vec = NULL;

	BGN_TRY_CATCH

	/** send signal */
	if (comm_rank_ == 0) {
		int sig = sig_addRow;
		MPI_Bcast(&sig, 1, MPI_INT, 0, comm_);
		DSPdebugMessage("Rank 0 sent signal %d.\n", sig);
		_size = vec->getNumElements();
		_indices = const_cast<int*>(vec->getIndices());
		_elements = const_cast<double*>(vec->getElements());
		_lb = lb;
		_ub = ub;
		_vec = const_cast<CoinPackedVector*>(vec);
	}

	/** synchronize information */
	MPI_Bcast(&_size, 1, MPI_INT, 0, comm_);
	DSPdebugMessage("Rank %d broadcasted size [%d].\n", comm_rank_, _size);
	if (comm_rank_ > 0) {
		_indices = new int [_size];
		_elements = new double [_size];
	}
	MPI_Bcast(_indices, _size, MPI_INT, 0, comm_);
	MPI_Bcast(_elements, _size, MPI_DOUBLE, 0, comm_);
	MPI_Bcast(&_lb, 1, MPI_DOUBLE, 0, comm_);
	MPI_Bcast(&_ub, 1, MPI_DOUBLE, 0, comm_);
	if (comm_rank_ > 0) {
		_vec = new CoinPackedVector(_size, _indices, _elements);
	}
	DwWorker::addRow(_vec, _lb, _ub);

	if (comm_rank_ == 0) {
		_indices = NULL;
		_elements = NULL;
		_vec = NULL;
	}

	END_TRY_CATCH(FREE_MEMORY)
	FREE_MEMORY
#undef FREE_MEMORY
}
