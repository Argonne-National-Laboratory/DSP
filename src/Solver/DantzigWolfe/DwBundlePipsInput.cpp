/*
 * DwBundlePipsInput.cpp
 *  Created on: July 16, 2018
 *      Author: Kibaek Kim
 */

//#define DSP_DEBUG
#include <numeric>
#include "Solver/DantzigWolfe/DwBundlePipsInput.h"
#include "Utility/DspMessage.h"

//#define SCALE_TAU

DwBundlePipsInput::DwBundlePipsInput(int nscen, int nvars1) : 
nscen_(nscen), 
nvars1_(nvars1),
ncons1_(0),
ncons2_(1),
tau_(1.0) {
	assert(nscen_ > 0);
	assert(nvars1_ > 0);
	bundle_x_.resize(nscen_);
	bundle_obj_.resize(nscen_);
}

int DwBundlePipsInput::nSecondStageVars(int scen) {
	return bundle_x_[scen].size();
}

std::vector<double> DwBundlePipsInput::getSecondStageColLB(int scen) {
	return std::vector<double>(nSecondStageVars(scen), 0.0);
}

std::vector<double> DwBundlePipsInput::getSecondStageColUB(int scen) {
	return std::vector<double>(nSecondStageVars(scen), COIN_DBL_MAX);
}

std::vector<double> DwBundlePipsInput::getSecondStageObj(int scen) {
	std::vector<double> obj; 
	obj.reserve(bundle_obj_[scen].size());
	double c;

	for (size_t j = 0; j < bundle_obj_[scen].size(); ++j) {
		c = bundle_obj_[scen][j];
		// if (scen == 0) printf("bundle_obj_[%d] = %e\n", j, bundle_obj_[j]);
		for (int k = 0; k < bundle_x_[scen][j].getNumElements(); ++k) {
			// if (scen == 0) 
			// 	printf("prox_center_[%d][%d] = %e, bundle_x_[%d].getElements()[%d] = %e\n", 
			// 		scen, bundle_x_[j].getIndices()[k], prox_center_[scen][bundle_x_[j].getIndices()[k]], 
			// 		j, k, bundle_x_[j].getElements()[k]);
			assert(bundle_x_[scen][j].getIndices()[k] < prox_center_[scen].size());
			c -= prox_center_[scen][bundle_x_[scen][j].getIndices()[k]] * bundle_x_[scen][j].getElements()[k];
		}
		// assert(prox_center_[scen].size() >= bundle_x_[j].getMaxIndex());
		// c -= bundle_x_[j].dotProduct(&prox_center_[scen][0]);
#ifdef SCALE_TAU
		c *= tau_;
#endif
		obj.push_back(c);
	}
#ifdef DSP_DEBUG
	printf("getSecondStageObj(%d) [size %u]: \n", scen, obj.size());
	DspMessage::printArray(obj.size(), &obj[0]);
#endif

	return obj;
}

std::vector<std::string> DwBundlePipsInput::getSecondStageColNames(int scen) {
	return std::vector<std::string>(nSecondStageVars(scen), "z");
}

CoinPackedMatrix DwBundlePipsInput::getFirstStageConstraints() {
	std::vector<CoinBigIndex> starts(nvars1_+1, 0);
	return CoinPackedMatrix(true, 0, nvars1_, 0, 0, 0, &starts[0], 0);
}

CoinPackedMatrix DwBundlePipsInput::getLinkingConstraints(int scen) {
	std::vector<CoinBigIndex> starts(nvars1_+1, 0);
	return CoinPackedMatrix(true, 1, nvars1_, 0, 0, 0, &starts[0], 0);
}

CoinPackedMatrix DwBundlePipsInput::getSecondStageConstraints(int scen) {
	CoinPackedMatrix W(true, 0, 0);
	W.setDimensions(ncons2_, 0);

	/** add columns */
	int nvars2 = nSecondStageVars(scen);
	int indices[] = {0};
	double elements[] = {1.0};
	for (size_t j = 0; j < nvars2; ++j)
		W.appendCol(1, indices, elements);
#ifdef DSP_DEBUG
	printf("getSecondStageConstraints(%d): \n", scen);
	DspMessage::printArray(W.getNumElements(), W.getIndices(), W.getElements());
#endif
	assert(W.verifyMtx(0) == 0);
	assert(W.getNumRows() == ncons2_);
	assert(W.getNumCols() == nvars2);

	return W;
}

CoinPackedMatrix DwBundlePipsInput::getFirstStageHessian() {
	CoinPackedMatrix H(true, 0, 0);
	H.setDimensions(nvars1_, 0);
	for (int i = 0; i < nvars1_; ++i) {
#ifdef SCALE_TAU
		double hesse = 1.0 * nscen_;
#else
		double hesse = 1.0 * nscen_ / tau_;
#endif
		H.appendCol(1, &i, &hesse);
	}
#ifdef DSP_DEBUG
	printf("getFirstStageHessian(): \n");
	DspMessage::printArray(H.getNumElements(), H.getIndices(), H.getElements());
#endif
	assert(H.verifyMtx(0) == 0);
	assert(H.getNumRows() == nvars1_);
	assert(H.getNumCols() == nvars1_);

	return H;
}

CoinPackedMatrix DwBundlePipsInput::getSecondStageCrossHessian(int scen) {
	CoinPackedMatrix H(false, 0, 0);
	H.setDimensions(0, nvars1_);

	/** add columns generated */
	std::vector<int> indices(nvars1_, 0);
	std::vector<double> values(nvars1_, 0.0);
	std::vector<double> xy;
	for (size_t s = 0; s < bundle_x_[scen].size(); ++s) {
		indices.clear(); values.clear();
		for (int j = 0; j < bundle_x_[scen][s].getNumElements(); ++j)
			if (fabs(bundle_x_[scen][s].getElements()[j]) > 1.0e-12) {
				indices.push_back(bundle_x_[scen][s].getIndices()[j]);
#ifdef SCALE_TAU
				values.push_back(-(bundle_x_[scen][s].getElements()[j]));
#else
				values.push_back(-(bundle_x_[scen][s].getElements()[j]) / tau_);
#endif
			}
		H.appendRow(indices.size(), &indices[0], &values[0]);
	}

	/** switch row-wise to column-wise for PIPS to understand */
	H.reverseOrdering();

#ifdef DSP_DEBUG
	printf("getSecondStageCrossHessian(%d): \n", scen);
	DspMessage::printArray(H.getNumElements(), H.getIndices(), H.getElements());
#endif
	assert(H.verifyMtx(0) == 0);
	assert(H.getNumRows() == nSecondStageVars(scen));
	assert(H.getNumCols() == nvars1_);

	return H;
}

CoinPackedMatrix DwBundlePipsInput::getSecondStageHessian(int scen) {
	/** number of second-stage variables */
	int nvars2 = nSecondStageVars(scen);

	CoinPackedMatrix H(true, 0, 0);
	H.setDimensions(nvars2, 0);

	/** add columns generated */
	std::vector<int> indices(nvars2, 0);
	std::vector<double> values(nvars2, 0.0);
	std::vector<double> dense_x(nvars1_, 0.0);
	for (size_t s1 = 0; s1 < bundle_x_[scen].size(); ++s1) {
		/** make it a dense vector */
		dense_x.assign(nvars1_, 0.0);
		for (int j = 0; j < bundle_x_[scen][s1].getNumElements(); ++j) {
			assert(bundle_x_[scen][s1].getIndices()[j] < nvars1_);
			dense_x[bundle_x_[scen][s1].getIndices()[j]] = bundle_x_[scen][s1].getElements()[j];
			// if (scen == 0) {
			// 	printf("bundle_x_[%d]: dense_x[%d] = %e\n", s1, bundle_x_[s1].getIndices()[j], dense_x[bundle_x_[s1].getIndices()[j]]);
			// }
		}
		/** compute x_i^T x_j / tau_ */
		indices.clear(); values.clear();
		for (size_t s2 = s1; s2 < bundle_x_[scen].size(); ++s2) {
			double xx = 0.0;
			for (int j = 0; j < bundle_x_[scen][s2].getNumElements(); ++j) {
				assert(bundle_x_[scen][s2].getIndices()[j] < nvars1_);
				// if (scen == 0) {
				// 	printf("%d/%d: x = %e, dense_x[%d] = %e, bundle_x_[%d] = %e\n", j, bundle_x_[s2].getNumElements(),
				// 		xx, bundle_x_[s2].getIndices()[j], dense_x[bundle_x_[s2].getIndices()[j]], 
				// 		bundle_x_[s2].getIndices()[j], bundle_x_[s2].getElements()[j]);
				// }
				xx += dense_x[bundle_x_[scen][s2].getIndices()[j]] * bundle_x_[scen][s2].getElements()[j];
				// if (scen == 0) printf("xx = %e\n", xx);
			}
			if (fabs(xx) > 1.0e-12) {
				indices.push_back(s2);
				// if (pos > H.getNumCols()) xx *= 2.0;
#ifdef SCALE_TAU						
				values.push_back(xx);
#else
				values.push_back(xx / tau_);
#endif						
			}
		}
		/** add the computed column */
		H.appendCol(indices.size(), &indices[0], &values[0]);
	}
#ifdef DSP_DEBUG
	printf("getSecondStageHessian(%d): \n", scen);
	DspMessage::printArray(H.getNumElements(), H.getIndices(), H.getElements());
#endif
	assert(H.verifyMtx(0) == 0);
	assert(H.getNumRows() == nvars2);
	assert(H.getNumCols() == nvars2);

	return H;
}
