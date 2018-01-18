/*
 * PipsInterface.cpp
 *  Created on: Dec 12, 2017
 *      Author: Kibaek Kim
 */

#include "Solver/DantzigWolfe/PipsInterface.h"
#include "Utility/DspMessage.h"
#include "PIPSIpmInterface.h"
//#include "sFactoryAug.h"
#include "sFactoryAugSchurLeaf.h"
#include "MehrotraStochSolver.h"

/** display option for PIPS */
int gOoqpPrintLevel = 0;

PipsInterface::PipsInterface(int nscen, int ncols) {
	input_ = new PipsInput(nscen, ncols);
	nrows_ = input_->nscen_ * input_->nvars1_;
	solution_.resize(input_->nscen_ + input_->nvars1_);
}

/** This function is call only by workers. */
void PipsInterface::addRow(const CoinPackedVector& v, double lb, double ub) {
	std::vector<int> cinds(v.getNumElements());
	std::vector<double> cvals(v.getNumElements());

	/** construct second-stage data */
	int sind = v.getIndices()[0];
	input_->matT_[sind].appendRow(0, NULL, NULL);
	cinds[0] = 0;
	cvals[0] = v.getElements()[0];
	for (int j = 1; j < v.getNumElements(); ++j) {
		cinds[j] = v.getIndices()[j] - input_->nscen_  + 1;
		cvals[j] = v.getElements()[j];
		//printf("cinds %d [%d]\n", j, cinds[j]);
	}
	input_->matW_[sind].appendRow(v.getNumElements(), &cinds[0], &cvals[0]);
	input_->rlbd2_[sind].push_back(lb);
	input_->rubd2_[sind].push_back(ub);
	input_->rname2_[sind].push_back("s" + std::to_string(sind) + "r" + std::to_string(input_->ncons2_[sind]));
	input_->ncons2_[sind]++;
	nrows_++;
	//input_->matT_[sind].verifyMtx(4);
	//input_->matW_[sind].verifyMtx(4);
}

/** This function is call only by workers. */
void PipsInterface::clearMatricesTW() {
	int ndels;
	std::vector<int> delinds;
	for (int j = 0; j < input_->nscen_; ++j) {
		ndels = input_->ncons2_[j] - input_->nvars1_;
		delinds.resize(ndels);
		std::iota(delinds.begin(), delinds.end(), 0);
		input_->matT_[j].deleteRows(ndels, &delinds[0]);
		input_->matW_[j].deleteRows(ndels, &delinds[0]);
		input_->rlbd2_[j].erase(input_->rlbd2_[j].begin() + input_->nvars1_, input_->rlbd2_[j].end());
		input_->rubd2_[j].erase(input_->rubd2_[j].begin() + input_->nvars1_, input_->rubd2_[j].end());
		input_->rname2_[j].erase(input_->rname2_[j].begin() + input_->nvars1_, input_->rname2_[j].end());
		input_->ncons2_[j] -= ndels;
		nrows_ -= ndels;
	}
}

void PipsInterface::retrieveBundles(std::vector<int>& rstart, std::vector<CoinPackedVector*>& vecs, std::vector<double>& rlbd, std::vector<double>& rubd, std::vector<int>& scen) {
	/** count the number of rows */
	int nrows = 0;
	for (int j = 0; j < input_->nscen_; ++j) {
		//printf("scen %d:\n", j);
		//input_->matW_[j].verifyMtx(4);
		nrows += input_->matW_[j].getNumRows() - rstart[j];
	}

	/** initialize data */
	vecs.clear(); rlbd.clear(); rubd.clear(); scen.clear();
	vecs.reserve(nrows); rlbd.reserve(nrows); rubd.reserve(nrows); scen.reserve(nrows);

	for (int j = 0; j < input_->nscen_; ++j) {
		const int* lens = input_->matW_[j].getVectorLengths();
		const CoinBigIndex* starts = input_->matW_[j].getVectorStarts();
		const int* indices = input_->matW_[j].getIndices();
		const double* elements = input_->matW_[j].getElements();
		for (int i = rstart[j]; i < input_->matW_[j].getNumRows(); ++i) {
			vecs.push_back(new CoinPackedVector(lens[i], indices + starts[i], elements + starts[i]));
			rlbd.push_back(input_->rlbd2_[j][i]);
			rubd.push_back(input_->rubd2_[j][i]);
			scen.push_back(j);
		}
		rstart[j] = input_->matW_[j].getNumRows(); 
	}
}

int PipsInterface::solve(double weight) {
	/** construct hessian */
	for (int i = 0; i < input_->nvars1_; ++i) 
		input_->hess1_.modifyCoefficient(i, i, weight);
	//input_->hess1_.verifyMtx(4);

	/** Create PIPS-IPM object */
	//printf("creating PIP-IPM object...\n");
	PIPSIpmInterface<sFactoryAugSchurLeaf, MehrotraStochSolver> solver(*input_);

	//printf("solving the master with PIPS-IPM...\n");
	solver.go();
	objval_ = solver.getObjective();

	std::vector<double> first_stage_solution = solver.getFirstStagePrimalColSolution();
	std::copy(first_stage_solution.begin(), first_stage_solution.end(), solution_.begin() + input_->nscen_);
	scenarios_.clear(); thetas_.clear();
	for (int j = 0; j < input_->nscen_; ++j) {
		std::vector<double> sol = solver.getSecondStagePrimalColSolution(j);
		if (sol.size() > 0) {
			//solution_[j] = sol[0];
			scenarios_.push_back(j);
			thetas_.push_back(sol[0]);
		}
	}

	return 0;
}

void PipsInterface::print() {
	printf("nscen_  %d\n", input_->nscen_);
	printf("nvars1_ %d\n", input_->nvars1_);
	printf("ncons1_ %d\n", input_->ncons1_);
	printf("nvars2_ %d\n", input_->nvars2_);
	printf("ncons2_");
	DspMessage::printArray(input_->nscen_, &(input_->ncons2_)[0]);
	//printf("clbd1_");
	//DspMessage::printArray(input_->nvars1_, &(input_->clbd1_)[0]);
	//printf("cubd1_");
	//DspMessage::printArray(input_->nvars1_, &(input_->cubd1_)[0]);
	printf("obj1_");
	DspMessage::printArray(input_->nvars1_, &(input_->obj1_)[0]);
	/*
	printf("cname1_");
	for (int j = 0; j < input_->nvars1_; ++j)
		printf(" %s", input_->cname1_[j].c_str());
	*/
	//printf("\nrlbd1_");
	//DspMessage::printArray(input_->ncons1_, &(input_->rlbd1_)[0]);
	//printf("rubd1_");
	//DspMessage::printArray(input_->ncons1_, &(input_->rubd1_)[0]);
	//printf("clbd2_");
	//DspMessage::printArray(input_->nvars2_, &(input_->clbd2_)[0]);
	//printf("cubd2_");
	//DspMessage::printArray(input_->nvars2_, &(input_->cubd2_)[0]);
	//printf("obj2_");
	/*
	for (int i = 0; i < input_->nscen_; ++i) {
		printf("scen %d:\n", i);
		printf("rlbd2_"); DspMessage::printArray(input_->ncons2_[i], &input_->rlbd2_[i][0]);
		printf("rubd2_"); DspMessage::printArray(input_->ncons2_[i], &input_->rubd2_[i][0]);
	}
	DspMessage::printArray(input_->nvars2_, &(input_->obj2_)[0]);
	for (int i = 0; i < input_->nscen_; ++i) {
		printf("\ncname2_[%d]", i);
		for (int j = 0; j < input_->nvars2_; ++j)
			printf(" %s", input_->cname2_[i][j].c_str());
	}
	*/
	printf("matA\n");
	input_->matA_.verifyMtx(4);
	for (int i = 0; i < input_->nscen_; ++i) {
		printf("scen %d:\n", i);
		input_->matW_[i].verifyMtx(4);
		input_->matT_[i].verifyMtx(4);
	}
}

