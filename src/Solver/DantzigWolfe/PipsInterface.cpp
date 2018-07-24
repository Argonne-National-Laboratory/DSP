/*
 * PipsInterface.cpp
 *  Created on: Dec 12, 2017
 *      Author: Kibaek Kim
 */

#define DSP_DEBUG

#include <numeric>

#include "Solver/DantzigWolfe/PipsInterface.h"
#include "Utility/DspMessage.h"
#include "PIPSIpmInterface.h"
#include "sFactoryAug.h"
// #include "sFactoryAugSchurLeaf.h"
#include "MehrotraStochSolver.h"

/** display option for PIPS */
int gOoqpPrintLevel = 0;

PipsInterface::PipsInterface(int nscen, int ncols) : 
penalty_(1.0), 
objval_(0.0), 
soln_w_(NULL) {
	input_ = new PipsInput(nscen, ncols);
	nrows_ = input_->nscen_ * input_->nvars1_;
	dualsol_.resize(input_->nscen_ + nrows_);
	bestdualsol_.resize(input_->nscen_ + nrows_);
}


PipsInterface::PipsInterface(const PipsInterface& rhs): 
nrows_(rhs.nrows_), 
penalty_(rhs.penalty_), 
scenarios_(rhs.scenarios_), 
objval_(rhs.objval_), 
bestdualsol_(rhs.bestdualsol_) {
	input_ = new PipsInput(*(rhs.input_));
	soln_w_ = new CoinPackedVector(rhs.soln_w_);
	for (int i = 0; i < rhs.soln_z_.size(); ++i)
		soln_z_.push_back(new CoinPackedVector(rhs.soln_z_[i]));
	for (int i = 0; i < rhs.soln_theta_.size(); ++i)
		soln_theta_.push_back(new CoinPackedVector(rhs.soln_theta_[i]));
}

PipsInterface::~PipsInterface() {
	delete input_;
	clearSolutions();
}

void PipsInterface::clearSolutions() {
	if (soln_w_ != NULL)
		delete soln_w_;
	for (int i = 0; i < soln_z_.size(); ++i) {
		if (soln_z_[i] != NULL)
			delete soln_z_[i];
	}
	soln_z_.clear();
	for (int i = 0; i < soln_theta_.size(); ++i) {
		if (soln_theta_[i] != NULL)
			delete soln_theta_[i];
	}
	soln_theta_.clear();
}

/** This function is called only by workers. */
void PipsInterface::addCol(int sind, const CoinPackedVector* x, double cx, std::vector<double>& bestdualsol) {
	int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	if (mype == 0) print();
	// DSPdebugMessage("Adding column...\n");

#ifdef DSP_DEBUG
	if (mype == 0) {
		printf("sind %d\n", sind);
		printf("cx: %e\nx:\n", cx);
		DspMessage::printArray(x);
	}
#endif

	int cinds[] = {0};
	double cvals[] = {1.0};
	input_->matW_[sind].appendCol(1, cinds, cvals);
	input_->clbd2_[sind].push_back(0.0);
	input_->cubd2_[sind].push_back(COIN_DBL_MAX);
	input_->cx2_[sind].push_back(cx);
	input_->x2_[sind].push_back(CoinPackedVector(*x));
	double obj2 = cx;
	for (int i = 0; i < x->getNumElements(); ++i)
		if (x->getIndices()[i] < nrows_) {
			int j = input_->nscen_ + sind * input_->nvars1_ + x->getIndices()[i] % input_->nvars1_;
			assert(j < bestdualsol.size());
			obj2 += bestdualsol[j] * x->getElements()[i];
		}
	input_->obj2_[sind].push_back(obj2);
	input_->cname2_[sind].push_back("s" + std::to_string(sind) + "r" + std::to_string(input_->nvars2_[sind]));
	input_->nvars2_[sind]++;

	/** add row to hessian matrix */
	std::vector<int> rowind;
	std::vector<double> rowelm;
	rowind.reserve(input_->nvars1_);
	rowelm.reserve(input_->nvars1_);
	double xy = 0.0, xval;
	for (int i = 0; i < x->getNumElements(); ++i)
		if (x->getIndices()[i] < nrows_) {
			rowind.push_back(x->getIndices()[i] % input_->nvars1_);
			xval = x->getElements()[i];
			rowelm.push_back(-xval/penalty_);
		}
	input_->hessxy_[sind].appendRow(rowind.size(), &rowind[0], &rowelm[0]);
	// input_->hessxy_[sind].appendCol(rowind.size(), &rowind[0], &rowelm[0]);
#ifdef DSP_DEBUG1
	if (mype == 0) {
		printf("new row to hessxy_:\n");
		DspMessage::printArray(rowind.size(), &rowind[0], &rowelm[0]);
	}
#endif

	rowind.clear();
	rowelm.clear();
	for (int i = 0; i < input_->x2_[sind].size(); ++i) {
		CoinPackedVector retval = input_->x2_[sind][i]*(*x);
		// for (int k = 0; k < input_->x2_[sind][i].getNumElements(); ++k)
		// 	printf("%d [%e] ", input_->x2_[sind][i].getIndices()[k], input_->x2_[sind][i].getElements()[k]);
		// printf("\n");
		// for (int k = 0; k < x->getNumElements(); ++k)
		// 	printf("%d [%e] ", x->getIndices()[k], x->getElements()[k]);
		// printf("\n");
		// for (int k = 0; k < retval.getNumElements(); ++k)
		// 	printf("%d [%e] ", retval.getIndices()[k], retval.getElements()[k]);
		// printf("\n");
		xy = 0.0;
		for(int j = 0; j < retval.getNumElements(); ++j){
			if (retval.getIndices()[j] < nrows_)
				xy += retval.getElements()[j];
		}
		if (fabs(xy) > 1.0e-8) {
			rowind.push_back(i);
			rowelm.push_back(xy/penalty_);
		}
	}
	input_->hessyy_[sind].setDimensions(input_->hessyy_[sind].getNumRows(), input_->x2_[sind].size());
	input_->hessyy_[sind].appendRow(rowind.size(), &rowind[0], &rowelm[0]);
#ifdef DSP_DEBUG1
	if (mype == 0) {
		printf("new row to hessyy_:\n");
		DspMessage::printArray(rowind.size(), &rowind[0], &rowelm[0]);
	}
#endif
	// input_->hessyy_[sind].setDimensions(input_->hessyy_[sind].getNumRows(), input_->x2_[sind].size());
	// input_->hessyy_[sind].appendRow(rowind.size(), &rowind[0], &rowelm[0]);
	// input_->hessxy_[sind].verifyMtx(4);

	// input_->hessxy_[sind].verifyMtx(4);
	// for (int i = 0; i < input_->hessxy_[sind].getNumRows(); ++i) {
	// 	int start = input_->hessxy_[sind].getVectorFirst(i);
	// 	int len = input_->hessxy_[sind].getVectorSize(i);
	// 	printf("row %d:\n", i);
	// 	printf("  indices:\n");
	// 	DspMessage::printArray(len, input_->hessxy_[sind].getIndices() + start);
	// 	printf("  elements:\n");
	// 	DspMessage::printArray(len, input_->hessxy_[sind].getElements() + start);
	// }

}

void PipsInterface::updateCenter(double penalty, std::vector<double>& bestdualsol) {
	DSPdebugMessage("Update center...\n");
	/** update linear coefficient */
	double obj2;
	for (int s = 0; s < input_->nscen_; ++s) {
		for (int k = 0; k < input_->cx2_[s].size(); ++k) {
			obj2 = input_->cx2_[s][k];
			for (int i = 0; i < input_->x2_[s][k].getNumElements(); ++i)
				if (input_->x2_[s][k].getIndices()[i] < input_->nvars1_ * input_->nscen_) {
					int j = input_->nscen_ + s * input_->nvars1_ + input_->x2_[s][k].getIndices()[i] % input_->nvars1_;
					if (j < bestdualsol.size())
						obj2 += bestdualsol[j] * input_->x2_[s][k].getElements()[i];
					else
						printf("Invalid indices (%d) for bestdualsol\n", j);
				}
			input_->obj2_[s][k] = obj2;
			DSPdebugMessage("obj2s_[%d][%d] = %e\n", s, k, input_->obj2_[s][k]);
		}
	}
	bestdualsol_ = bestdualsol;

	DSPdebugMessage("diff on penalty %e\n", fabs(penalty-penalty_));
	if (fabs(penalty-penalty_) > 1.0e-8) {
		/** adjust hessian */
		for (int i = 0; i < input_->hessxx_.getNumElements(); ++i) {
			//DSPdebugMessage("updating HessianXX %d [%e]\n", input_->hessxx_.getIndices()[i], input_->hessxx_.getElements()[i]);
			input_->hessxx_.getMutableElements()[i] *=  penalty_ / penalty;
		}
		for (int s = 0; s < input_->nscen_; ++s) {
			for (int i = 0; i < input_->hessxy_[s].getNumElements(); ++i) {
				DSPdebugMessage("updating HessianXY %d [%e]\n", input_->hessxy_[s].getIndices()[i], input_->hessxy_[s].getElements()[i]);
				input_->hessxy_[s].getMutableElements()[i] *=  penalty_ / penalty;
			}
			for (int i = 0; i < input_->hessyy_[s].getNumElements(); ++i) {
				//DSPdebugMessage("updating HessianXY %d [%e]\n", input_->hessyy_[s].getIndices()[i], input_->hessyy_[s].getElements()[i]);
				input_->hessyy_[s].getMutableElements()[i] *=  penalty_ / penalty;
			}
		}
		/** update penalty parameter */
		penalty_ = penalty;
	}
}

/** This function is call only by workers. */
void PipsInterface::clearMatricesW() {
	int ndels;
	std::vector<int> delinds;
	for (int j = 0; j < input_->nscen_; ++j) {
		ndels = input_->matW_[j].getNumCols();
		delinds.resize(ndels);
		std::iota(delinds.begin(), delinds.end(), 0);
		input_->matW_[j].deleteCols(ndels, &delinds[0]);
		input_->clbd2_[j].clear();
		input_->cubd2_[j].clear();
		input_->cx2_[j].clear();
		input_->x2_[j].clear();
		input_->obj2_[j].clear();
		input_->cname2_[j].clear();
		input_->nvars2_[j] = 0;
		ndels = input_->hessxy_[j].getNumRows();
		delinds.resize(ndels);
		std::iota(delinds.begin(), delinds.end(), 0);
		input_->hessxy_[j].deleteRows(ndels, &delinds[0]);
		ndels = input_->hessyy_[j].getNumRows();
		delinds.resize(ndels);
		std::iota(delinds.begin(), delinds.end(), 0);
		input_->hessyy_[j].deleteRows(ndels, &delinds[0]);
		input_->hessyy_[j].deleteCols(ndels, &delinds[0]);
	}
}

void PipsInterface::collectColumnData(
		double& penalty,
		std::vector<int>& cstart,
		std::vector<int>& scen,
		std::vector<CoinPackedVector*>& x,
		std::vector<double>& cx) {
	penalty = penalty_;
	/** count the number of columns to retrieve */
	int ncols = 0;
	std::vector<int> ncolsPerScen(input_->nscen_, 0);
	for (int j = 0; j < input_->nscen_; ++j) {
		ncolsPerScen[j] = input_->obj2_[j].size() - cstart[j];
		ncols += ncolsPerScen[j];
	}

	/** initialize data */
	scen.clear(); x.clear(); cx.clear();
	scen.reserve(ncols); x.reserve(ncols); cx.reserve(ncols);

	for (int j = 0; j < input_->nscen_; ++j) {
		for (int i = 0; i < ncolsPerScen[j]; ++i) {
			/** scenarion index */
			scen.push_back(j);
			/** objective coefficient for new column */
			cx.push_back(input_->cx2_[j][cstart[j] + i]);
			/** subproblem solution x */
			x.push_back(new CoinPackedVector(input_->x2_[j][cstart[j] + i]));
		}

		/** increment start column */
		cstart[j] += ncolsPerScen[j];
	}
}

int PipsInterface::solve(double weight) {

	// DEBUG
	int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	if (mype == 0) print();

	/** Create PIPS-IPM object */
	if (mype == 0) printf("creating PIPS-IPM object...\n");
	PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> solver(*input_);
	// PIPSIpmInterface<sFactoryAugSchurLeaf, MehrotraStochSolver> solver(*input_);

	if (mype == 0) printf("solving the master with PIPS-IPM...\n");
	solver.go();
	objval_ = solver.getObjective();
	if (mype == 0) printf("done with objective %e\n", objval_);

	/** clear solution storage */
	clearSolutions();

	/** store first-stage variable solution */
	std::vector<double> soln_w = solver.getFirstStagePrimalColSolution();
	assert(soln_w.size() == input_->nvars1_);
#ifdef DSP_DEBUG
	if (mype == 0) {
		printf("soln w:\n");
		DspMessage::printArray(soln_w.size(), &soln_w[0]);
	}
#endif
	soln_w_ = new CoinPackedVector;
	for (int i = 0; i < soln_w.size(); ++i)
		if (fabs(soln_w[i]) > 1.0e-8)
			soln_w_->insert(i, soln_w[i]);

	std::vector<double> lambda(input_->nvars1_, 0.0);
	std::vector<double> beta(input_->nvars1_, 0.0);
	scenarios_.clear();

	int comm_size; MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	for (int r = 0; r < comm_size; ++r) {
		if (r == mype) {

	for (int s = 0; s < input_->nscen_; ++s) {

		std::vector<double> soln_z = solver.getSecondStagePrimalColSolution(s);
		if (soln_z.size() > 0) {
			assert(soln_z.size() == input_->x2_[s].size());
#ifdef DSP_DEBUG
			printf("soln z[%d]:\n", s);
			DspMessage::printArray(soln_z.size(), &soln_z[0]);
#endif
			/**  calculating lambda */
				// lambda[i] = bestdualsol_[input_->nscen_ + input_->nvars1_ * s + i] - soln_w[i] / penalty_;
			for (int i = 0; i < input_->nvars1_; ++i)
				beta[i] = -soln_w[i] / penalty_;
			for (int i = 0; i < input_->x2_[s].size(); ++i)
				for (int j = 0; j < input_->x2_[s][i].getNumElements(); ++j) {
					if (input_->x2_[s][i].getIndices()[j] < input_->nvars1_ * input_->nscen_) {
						beta[input_->x2_[s][i].getIndices()[j] % input_->nvars1_] += soln_z[i] * input_->x2_[s][i].getElements()[j] / penalty_;
#ifdef DSP_DEBUG
						// if (input_->x2_[s][i].getIndices()[j] % input_->nvars1_ == 0)
							printf("soln_z %e x2[%d] %e\n", soln_z[i], input_->x2_[s][i].getIndices()[j], input_->x2_[s][i].getElements()[j]);
#endif
					}
				}
#ifdef DSP_DEBUG
			printf("beta[%d]:\n", s);
			DspMessage::printArray(beta.size(), &beta[0]);
#endif
			for (int i = 0; i < input_->nvars1_; ++i) {
				lambda[i] = beta[i] + bestdualsol_[input_->nscen_ + input_->nvars1_ * s + i];
				//DSPdebugMessage("0.5 * tau * beta^2[%d] = %e\n", i, beta[i]*beta[i] * penalty_ / 2);
			}

			CoinPackedVector* sparse_lambda = new CoinPackedVector;
			for (int i = 0; i < lambda.size(); ++i)
				if (fabs(lambda[i]) > 1.0e-8)
					sparse_lambda->insert(i, lambda[i]);
			soln_lambda_.push_back(sparse_lambda);
			sparse_lambda = NULL;

			/**  calculating theta */
			std::vector<double> soln_theta = solver.getSecondStageDualRowSolution(s);
#ifdef DSP_DEBUG
			printf("theta[%d]:\n", s);
			DspMessage::printArray(soln_theta.size(), &soln_theta[0]);
#endif
			CoinPackedVector* sparse_theta = new CoinPackedVector;
			for (int i = 0; i < soln_theta.size(); ++i)
				if (fabs(soln_theta[i]) > 1.0e-8)
					sparse_theta->insert(i, soln_theta[i]);
			soln_theta_.push_back(sparse_theta);
			sparse_theta = NULL;

			/** store scenario index */
			scenarios_.push_back(s);
		}
		// if (solz.size() == input_->x2_[s].size()) {
		// 	for (int k = 0; k < input_->x2_[s].size(); ++k)
		// 		for (int i = 0; i < input_->x2_[s][k].getNumElements(); ++i) {
		// 			int j = s * input_->nvars1_ + input_->x2_[s][k].getIndices()[i] % input_->nvars1_;
		// 			lambda_[j] += solz[k] * input_->x2_[s][k].getElements()[i] / penalty_;
		// 		}
		// } else {
		// 	printf("solz and x2_[%d] do not match dimentsions.\n", s);
		// }
	}
	

		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	return 0;
}

void PipsInterface::print() {
	printf("nscen_  %d\n", input_->nscen_);
	printf("nvars1_ %d\n", input_->nvars1_);
	printf("ncons1_ %d\n", input_->ncons1_);
	printf("nvars2_\n");
	DspMessage::printArray(input_->nvars2_.size(), &(input_->nvars2_)[0]);
	printf("ncons2_ %d\n", input_->ncons2_);
	printf("obj2_\n");
	for (int s = 0; s < input_->nscen_; ++s) {
		printf("  scen %d:\n", s);
		DspMessage::printArray(input_->nvars2_[s], &(input_->obj2_[s])[0]);
	}
	// printf("matA_\n");
	// input_->matA_.verifyMtx(4);
	printf("hessxx_\n");
	input_->hessxx_.verifyMtx(4);
	for (int j = 0; j < input_->hessxx_.getNumCols(); ++j) {
		printf(" col %d\n", j);
		DspMessage::printArray(input_->hessxx_.getVectorSize(j), 
			input_->hessxx_.getIndices() + input_->hessxx_.getVectorFirst(j),
			input_->hessxx_.getElements() + input_->hessxx_.getVectorFirst(j));
	}
	for (int i = 0; i < input_->nscen_; ++i) {
		printf("scen %d:\n", i);
		// printf("\n  matW_:\n", i);
		// input_->matW_[i].verifyMtx(4);
		// DspMessage::printArray(input_->matW_[i].getNumCols(), input_->matW_[i].getVectorStarts());
		// DspMessage::printArray(input_->matW_[i].getNumCols(), input_->matW_[i].getVectorLengths());
		// DspMessage::printArray(input_->matW_[i].getNumElements(), input_->matW_[i].getIndices());
		// DspMessage::printArray(input_->matW_[i].getNumElements(), input_->matW_[i].getElements());
		printf("\n  hessxy_:\n", i);
		input_->hessxy_[i].verifyMtx(4);
		for (int j = 0; j < input_->hessxy_[i].getNumCols(); ++j) {
			printf(" col %d\n", j);
			DspMessage::printArray(input_->hessxy_[i].getVectorSize(j), 
				input_->hessxy_[i].getIndices() + input_->hessxy_[i].getVectorFirst(j),
				input_->hessxy_[i].getElements() + input_->hessxy_[i].getVectorFirst(j));
		}
		printf("\n  hessyy_:\n", i);
		input_->hessyy_[i].verifyMtx(4);
		for (int j = 0; j < input_->hessyy_[i].getNumCols(); ++j) {
			printf(" col %d\n", j);
			DspMessage::printArray(input_->hessyy_[i].getVectorSize(j), 
				input_->hessyy_[i].getIndices() + input_->hessyy_[i].getVectorFirst(j),
				input_->hessyy_[i].getElements() + input_->hessyy_[i].getVectorFirst(j));
		}
	}
}

