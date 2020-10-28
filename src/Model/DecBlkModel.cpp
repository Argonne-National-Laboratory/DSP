/*
 * DecBlkModel.cpp
 *
 *  Created on: Aug 31, 2016
 *      Author: kibaekkim
 */

//#define DSP_DEBUG

#include <vector>
/** Coin */
#include "CoinPackedVector.hpp"
/** Dsp */
#include "Utility/DspMessage.h"
#include "Model/DecBlkModel.h"

DecBlkModel::DecBlkModel() {
	blk_ = new BlkModel();
}

DecBlkModel::DecBlkModel(const DecBlkModel& rhs) {
	BlkModel* blk = const_cast<DecBlkModel&>(rhs).blkPtr();
	blk_ = new BlkModel(*blk);
}

DecBlkModel::~DecBlkModel() {
	FREE_PTR(blk_);
}

double DecBlkModel::evalLhsCouplingRow(int row, double** solutions) {
	double val = 0.0;
	for (int i = 0; i < blk_->getNumBlocks() - 1; ++i)
		val += evalLhsCouplingRowSubprob(row, i, solutions[i]);
	return val;
}

double DecBlkModel::evalLhsCouplingRowSubprob(
		int row,
		int subprob,
		double* subprobSolution) {
	/** retrieve block */
	DetBlock* block = blk_->block(subprob+1);
	assert(block != NULL);
	assert(row >= 0 && row < block->getNumCouplingRows());

	/** retrieve the master matrix */
	const CoinPackedMatrix* mat = blk_->block(0)->getConstraintMatrix();

	/** calculate lhs */
	double val = 0.0;
	for (int j = 0; j < block->getNumCouplingCols(); ++j) {
		/** get column index w.r.t. the master */
		int col = block->getCouplingCols()[j];
		val += mat->getCoefficient(row,col) * subprobSolution[col];

	}
	return val;
}

char DecBlkModel::getSenseCouplingRow(int row) {
	/** retrieve block */
	DetBlock* block = blk_->block(0);
	/** row index w.r.t. the master */
	int i = block->getCouplingRows()[row];
	DSPdebugMessage("row index %d -> %d\n", row, i);
	char sense = 'R';
	if (block->getRowLower()[i] < -1e+20)
		sense = 'L';
	else if (block->getRowUpper()[i] > +1e+20)
		sense = 'G';
	else if (block->getRowUpper()[i] - block->getRowLower()[i] < 1e-10)
		sense = 'E';
	else
		fprintf(stderr, "Coupling row cannot be described as a ranged constraint.\n");
	return sense;
}

double DecBlkModel::getRhsCouplingRow(int row) {
	/** retrieve block */
	DetBlock* block = blk_->block(0);
	/** row index w.r.t. the master */
	int i = block->getCouplingRows()[row];
	/** retrieve the master matrix */
	double rhs = 0.0;
	if (block->getRowLower()[i] < -1e+20)
		rhs = block->getRowUpper()[i];
	else if (block->getRowUpper()[i] > +1e+20)
		rhs = block->getRowLower()[i];
	else if (block->getRowUpper()[i] - block->getRowLower()[i] < 1e-10)
		rhs = block->getRowUpper()[i];
	else
		fprintf(stderr, "Coupling row cannot be described as a ranged constraint.\n");
	return rhs;
}

DSP_RTN_CODE DecBlkModel::decompose(
		int size,                /**< [in] size of subproblem subset */
		int * subprobs,          /**< [in] subset of subproblems */
		int naux,                /**< [in] number of auxiliary columns */
		double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) {

	BGN_TRY_CATCH

	/** the master problem part */
	DetBlock* master = blk_->block(0);
	/** number of rows and columns */
	int nrows = master->getNumRows();
	int ncols = master->getNumCols();
	/** constraint matrix */
	mat = new CoinPackedMatrix(*(master->getConstraintMatrix()));
	/** make sure the matrix is row-wise */
	if (mat->isColOrdered()) mat->reverseOrdering();

	/** the subproblem part */
	int coloffset = 0;
	for (int i = 0; i < size; ++i) {
		/** retrieve sub-block */
		DetBlock* sub = blk_->block(subprobs[i]+1);
		/** number of rows and columns */
		nrows += sub->getNumRows();
		ncols += sub->getNumCols();
		/** constraint matrix */
		const CoinPackedMatrix* submat = sub->getConstraintMatrix();
		const int* indices = submat->getIndices();
		const int* start = submat->getVectorStarts();
		const double* values = submat->getElements();
		/** add rows with column indices adjusted */
		for (int k = 0; k < submat->getNumRows(); ++k) {
			CoinPackedVector row;
			for (int j = 0; j < submat->getVectorSize(k); ++j) {
				if (indices[start[k]+j] < mat->getNumCols())
					row.insert(indices[start[k]+j], values[start[k]+j]);
				else
					row.insert(indices[start[k]+j] + coloffset, values[start[k]+j]);
			}
			mat->appendRow(row);
		}
		coloffset += sub->getNumCols() - mat->getNumCols();
		sub = NULL;
	}

	/** auxiliary variables */
	ncols += naux;

	/** allocate memory for the outputs */
	clbd  = new double [ncols];
	cubd  = new double [ncols];
	ctype = new char [ncols];
	obj   = new double [ncols];
	rlbd  = new double [nrows];
	rubd  = new double [nrows];

	/** copy the others */
	CoinCopyN(master->getColLower(), master->getNumCols(), clbd);
	CoinCopyN(master->getColUpper(), master->getNumCols(), cubd);
	CoinCopyN(master->getCtype(), master->getNumCols(), ctype);
	CoinCopyN(master->getObj(), master->getNumCols(), obj);
	CoinCopyN(master->getRowLower(), master->getNumRows(), rlbd);
	CoinCopyN(master->getRowUpper(), master->getNumRows(), rubd);
	/** mark positions */
	int cpos = master->getNumCols();
	int rpos = master->getNumRows();
	/** for sub-blocks */
	for (int i = 0; i < size; ++i) {
		DetBlock* sub = blk_->block(subprobs[i]+1);
		CoinCopyN(sub->getColLower(), sub->getNumCols(), clbd + cpos);
		CoinCopyN(sub->getColUpper(), sub->getNumCols(), cubd + cpos);
		CoinCopyN(sub->getCtype(), sub->getNumCols(), ctype + cpos);
		CoinCopyN(sub->getObj(), sub->getNumCols(), obj + cpos);
		CoinCopyN(sub->getRowLower(), sub->getNumRows(), rlbd + rpos);
		CoinCopyN(sub->getRowUpper(), sub->getNumRows(), rubd + rpos);
		cpos += sub->getNumCols();
		rpos += sub->getNumRows();
		sub = NULL;
	}
	/** for auxiliary variables */
	CoinCopyN(clbd_aux, naux, clbd + cpos);
	CoinCopyN(cubd_aux, naux, cubd + cpos);
	CoinFillN(ctype + cpos, naux, 'C');
	CoinCopyN(obj_aux, naux, obj + cpos);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DecBlkModel::decompose(
		int size,                /**< [in] size of subproblem subset */
		int * subprobs,          /**< [in] subset of subproblems */
		int naux,                /**< [in] number of auxiliary columns */
		double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		CoinPackedMatrix *& qobj,/**< [out] quadratic objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) {

	BGN_TRY_CATCH

	/** the master problem part */
	DetBlock* master = blk_->block(0);
	/** number of rows and columns */
	int nrows = master->getNumRows();
	int ncols = master->getNumCols();
	/** constraint matrix */
	mat = new CoinPackedMatrix(*(master->getConstraintMatrix()));
	/** make sure the matrix is row-wise */
	if (mat->isColOrdered()) mat->reverseOrdering();

	/** master quadratic objective */
	CoinPackedMatrix *qmaster=NULL;
	qmaster = new CoinPackedMatrix(*master->getQuadraticObjectiveMatrix());
	if (qmaster->isColOrdered()) qmaster->reverseOrdering();

	vector<int> qrowIndices;
	vector<int> qcolIndices;
	vector<double> qelements;

	int rownum=0;
	int pos=0;

	if (qmaster !=NULL){
		const CoinBigIndex * qstarts=qmaster->getVectorStarts();
		for (int i = 0; i < qmaster->getNumRows(); ++i)
		{
			int qlength = qmaster->getVectorSize(i);
			DSPdebugMessage("length of the %d th vector = %d", i, qlength);
			for (int j = 0; j < qlength; ++j) {
					qrowIndices.push_back(rownum);
					qcolIndices.push_back(qmaster->getIndices()[qstarts[i]+j]);
					if (rownum!=qmaster->getIndices()[qstarts[i]+j]){
						qelements.push_back(2*qmaster->getElements()[qstarts[i]+j]);
					}
					else{
						qelements.push_back(qmaster->getElements()[qstarts[i]+j]);
					}

					pos++;
			}	
			rownum++;
		}
	}


	/** the subproblem part */
	int coloffset = 0;
	for (int i = 0; i < size; ++i) {
		/** retrieve sub-block */
		DetBlock* sub = blk_->block(subprobs[i]+1);
		/** number of rows and columns */
		nrows += sub->getNumRows();
		ncols += sub->getNumCols();
		/** constraint matrix */
		const CoinPackedMatrix* submat = sub->getConstraintMatrix();
		const int* indices = submat->getIndices();
		const int* start = submat->getVectorStarts();
		const double* values = submat->getElements();

		/** quadratic objective */
		const CoinPackedMatrix* subqobj = sub->getQuadraticObjectiveMatrix();

		/** add rows with column indices adjusted */
		for (int k = 0; k < submat->getNumRows(); ++k) {
			CoinPackedVector row;
			for (int j = 0; j < submat->getVectorSize(k); ++j) {
				if (indices[start[k]+j] < mat->getNumCols())
					row.insert(indices[start[k]+j], values[start[k]+j]);
				else
					row.insert(indices[start[k]+j] + coloffset, values[start[k]+j]);
			}
			mat->appendRow(row);
		}

		if (subqobj != NULL){
			int rownum=0;
			for (int k = 0; k < subqobj->getNumRows(); ++k){
				CoinBigIndex qstart = subqobj->getVectorStarts()[k];
				for (int j=0; j < subqobj->getVectorSize(k); j++){
					if (subqobj->getIndices()[qstart+j]<mat->getNumCols()){
						qcolIndices.push_back(subqobj->getIndices()[qstart+j]);
					}
					else{
						qcolIndices.push_back(subqobj->getIndices()[qstart+j]+coloffset);
					}
					if(rownum<mat->getNumCols()){
						qrowIndices.push_back(rownum);
					}
					else{
						qrowIndices.push_back(rownum+coloffset);
					}
					if (rownum!=subqobj->getIndices()[qstart+j]){
						qelements.push_back(2 * (subqobj->getElements()[qstart + j]));
					}
					else{
						qelements.push_back((subqobj->getElements()[qstart + j]));
					}
				}
				rownum++;
			}
		}
		coloffset += sub->getNumCols() - mat->getNumCols();
		sub = NULL;
	}
	qobj=new CoinPackedMatrix(false, &qrowIndices[0], &qcolIndices[0], &qelements[0], qrowIndices.size());
	
	/** auxiliary variables */
	ncols += naux;

	/** allocate memory for the outputs */
	clbd  = new double [ncols];
	cubd  = new double [ncols];
	ctype = new char [ncols];
	obj   = new double [ncols];
	rlbd  = new double [nrows];
	rubd  = new double [nrows];

	/** copy the others */
	CoinCopyN(master->getColLower(), master->getNumCols(), clbd);
	CoinCopyN(master->getColUpper(), master->getNumCols(), cubd);
	CoinCopyN(master->getCtype(), master->getNumCols(), ctype);
	CoinCopyN(master->getObj(), master->getNumCols(), obj);
	CoinCopyN(master->getRowLower(), master->getNumRows(), rlbd);
	CoinCopyN(master->getRowUpper(), master->getNumRows(), rubd);
	/** mark positions */
	int cpos = master->getNumCols();
	int rpos = master->getNumRows();
	/** for sub-blocks */
	for (int i = 0; i < size; ++i) {
		DetBlock* sub = blk_->block(subprobs[i]+1);
		CoinCopyN(sub->getColLower(), sub->getNumCols(), clbd + cpos);
		CoinCopyN(sub->getColUpper(), sub->getNumCols(), cubd + cpos);
		CoinCopyN(sub->getCtype(), sub->getNumCols(), ctype + cpos);
		CoinCopyN(sub->getObj(), sub->getNumCols(), obj + cpos);
		CoinCopyN(sub->getRowLower(), sub->getNumRows(), rlbd + rpos);
		CoinCopyN(sub->getRowUpper(), sub->getNumRows(), rubd + rpos);
		cpos += sub->getNumCols();
		rpos += sub->getNumRows();
		sub = NULL;
	}
	/** for auxiliary variables */
	CoinCopyN(clbd_aux, naux, clbd + cpos);
	CoinCopyN(cubd_aux, naux, cubd + cpos);
	CoinFillN(ctype + cpos, naux, 'C');
	CoinCopyN(obj_aux, naux, obj + cpos);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}


DSP_RTN_CODE DecBlkModel::copySubprob(
		int subprob,             /**< [in] subproblem index */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix [A^k B^k] */
		double *& clbd,          /**< [out] column lower bounds of y */
		double *& cubd,          /**< [out] column upper bounds of y */
		char   *& ctype,         /**< [out] column types of y */
		double *& obj,           /**< [out] objective coefficients for y */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) {
	BGN_TRY_CATCH

	/** retrieve sub-block */
	DetBlock* sub = blk_->block(subprob+1);

	/** copy matrix */
	mat = new CoinPackedMatrix(*(sub->getConstraintMatrix()));
	DSPdebug(mat->verifyMtx(4));

	/** allocate memory */
	clbd = new double [mat->getNumCols()];
	cubd = new double [mat->getNumCols()];
	ctype = new char [mat->getNumCols()];
	obj = new double [mat->getNumCols()];
	rlbd = new double [mat->getNumRows()];
	rubd = new double [mat->getNumRows()];
	DSPdebugMessage("ncols %d nrows %d\n", mat->getNumCols(), mat->getNumRows());

	/** copy data */
	CoinCopyN(sub->getColLower(), sub->getNumCols(), clbd);
	CoinCopyN(sub->getColUpper(), sub->getNumCols(), cubd);
	CoinCopyN(sub->getCtype(), sub->getNumCols(), ctype);
	CoinCopyN(sub->getObj(), sub->getNumCols(), obj);
	CoinCopyN(sub->getRowLower(), sub->getNumRows(), rlbd);
	CoinCopyN(sub->getRowUpper(), sub->getNumRows(), rubd);
#ifdef DSP_DEBUG1
	printf("clbd:\n"); DspMessage::printArray(sub->getNumCols(), sub->getColLower());
	printf("cubd:\n"); DspMessage::printArray(sub->getNumCols(), sub->getColUpper());
	printf("obj:\n");  DspMessage::printArray(sub->getNumCols(), sub->getObj());
	printf("rlbd:\n"); DspMessage::printArray(sub->getNumRows(), sub->getRowLower());
	printf("rubd:\n"); DspMessage::printArray(sub->getNumRows(), sub->getRowUpper());
#endif

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DecBlkModel::copyRecoProb(
		int scen,                     /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech, /**< [out] technology matrix (A matrix) */
		CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix (B matrix) */
		double *& clbd_reco,          /**< [out] column lower bounds of y */
		double *& cubd_reco,          /**< [out] column upper bounds of y */
		char   *& ctype_reco,         /**< [out] column types of y */
		double *& obj_reco,           /**< [out] objective coefficients for y */
		double *& rlbd_reco,          /**< [out] row lower bounds */
		double *& rubd_reco,          /**< [out] row upper bounds */
		bool adjust_probability       /**< not used */) {

	BGN_TRY_CATCH

	/** the master problem part */
	DetBlock* master = blk_->block(0);
	/** number of rows and columns */
	int nrows = master->getNumRows();
	int ncols = master->getNumCols();
	/** constraint matrix */
	CoinPackedMatrix* mat = new CoinPackedMatrix(*(master->getConstraintMatrix()));

	/** allocate memory */
	mat_tech = new CoinPackedMatrix(false, 0, 0);
	mat_reco = new CoinPackedMatrix(false, 0, 0);

	/** retrieve sub-block */
	DetBlock* sub = blk_->block(scen+1);
	/** number of rows and columns */
	nrows += sub->getNumRows();
	ncols += sub->getNumCols();
	/** constraint matrix */
	const CoinPackedMatrix* submat = sub->getConstraintMatrix();
	const int* indices = submat->getIndices();
	const int* start = submat->getVectorStarts();
	const double* values = submat->getElements();
	/** resize matrices */
	mat_tech->setDimensions(0, mat->getNumCols());
	mat_reco->setDimensions(0, submat->getNumCols());
	/** add rows with column indices adjusted */
	for (int k = 0; k < submat->getNumRows(); ++k) {
		CoinPackedVector row_tech, row_reco;
		for (int j = 0; j < submat->getVectorSize(k); ++j) {
			if (indices[start[k]+j] < mat->getNumCols())
				row_tech.insert(indices[start[k]+j], values[start[k]+j]);
			else
				row_reco.insert(indices[start[k]+j] - mat->getNumCols(), values[start[k]+j]);
		}
		mat_tech->appendRow(row_tech);
		mat_reco->appendRow(row_reco);
	}

	/** allocate memory for the outputs */
	clbd_reco  = new double [sub->getNumCols()];
	cubd_reco  = new double [sub->getNumCols()];
	ctype_reco = new char [sub->getNumCols()];
	obj_reco   = new double [sub->getNumCols()];
	rlbd_reco  = new double [sub->getNumRows()];
	rubd_reco  = new double [sub->getNumRows()];

	/** copy data */
	CoinCopyN(sub->getColLower(), sub->getNumCols(), clbd_reco);
	CoinCopyN(sub->getColUpper(), sub->getNumCols(), cubd_reco);
	CoinCopyN(sub->getCtype(), sub->getNumCols(), ctype_reco);
	CoinCopyN(sub->getObj(), sub->getNumCols(), obj_reco);
	CoinCopyN(sub->getRowLower(), sub->getNumRows(), rlbd_reco);
	CoinCopyN(sub->getRowUpper(), sub->getNumRows(), rubd_reco);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DecBlkModel::copyRecoProb(
		int scen,                     /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech, /**< [out] technology matrix (A matrix) */
		CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix (B matrix) */
		double *& clbd_reco,          /**< [out] column lower bounds of y */
		double *& cubd_reco,          /**< [out] column upper bounds of y */
		char   *& ctype_reco,         /**< [out] column types of y */
		double *& obj_reco,           /**< [out] objective coefficients for y */
		CoinPackedMatrix *& qobj_reco_coupling,/**< [out] coupling quadratric coefficients (y^2)*/
		CoinPackedMatrix *& qobj_reco_ncoupling, /**< [out] non-coupling quadratic coefficients (xy) */
		double *& rlbd_reco,          /**< [out] row lower bounds */
		double *& rubd_reco           /**< [out] row upper bounds */) {

	BGN_TRY_CATCH

	/** the master problem part */
	DetBlock* master = blk_->block(0);
	/** number of rows and columns */
	int nrows = master->getNumRows();
	int ncols = master->getNumCols();
	/** constraint matrix */
	CoinPackedMatrix* mat = new CoinPackedMatrix(*(master->getConstraintMatrix()));

	/** allocate memory */
	mat_tech = new CoinPackedMatrix(false, 0, 0);
	mat_reco = new CoinPackedMatrix(false, 0, 0);

	/** retrieve sub-block */
	DetBlock* sub = blk_->block(scen+1);
	/** number of rows and columns */
	nrows += sub->getNumRows();
	ncols += sub->getNumCols();
	/** constraint matrix */
	const CoinPackedMatrix* submat = sub->getConstraintMatrix();
	const int* indices = submat->getIndices();
	const int* start = submat->getVectorStarts();
	const double* values = submat->getElements();
	/** resize matrices */
	mat_tech->setDimensions(0, mat->getNumCols());
	mat_reco->setDimensions(0, submat->getNumCols());
	/** add rows with column indices adjusted */
	for (int k = 0; k < submat->getNumRows(); ++k) {
		CoinPackedVector row_tech, row_reco;
		for (int j = 0; j < submat->getVectorSize(k); ++j) {
			if (indices[start[k]+j] < mat->getNumCols())
				row_tech.insert(indices[start[k]+j], values[start[k]+j]);
			else
				row_reco.insert(indices[start[k]+j] - mat->getNumCols(), values[start[k]+j]);
		}
		mat_tech->appendRow(row_tech);
		mat_reco->appendRow(row_reco);
	}

	/** allocate memory for the outputs */
	clbd_reco  = new double [sub->getNumCols()];
	cubd_reco  = new double [sub->getNumCols()];
	ctype_reco = new char [sub->getNumCols()];
	obj_reco   = new double [sub->getNumCols()];
	rlbd_reco  = new double [sub->getNumRows()];
	rubd_reco  = new double [sub->getNumRows()];

	/** copy data */
	CoinCopyN(sub->getColLower(), sub->getNumCols(), clbd_reco);
	CoinCopyN(sub->getColUpper(), sub->getNumCols(), cubd_reco);
	CoinCopyN(sub->getCtype(), sub->getNumCols(), ctype_reco);
	CoinCopyN(sub->getObj(), sub->getNumCols(), obj_reco);
	CoinCopyN(sub->getRowLower(), sub->getNumRows(), rlbd_reco);
	CoinCopyN(sub->getRowUpper(), sub->getNumRows(), rubd_reco);
	//CoinPackedMatrix *qobj_reco;
	//qobj_reco->copyOf(*sub->getQuadraticObjectiveMatrix());
	// TO BE MODIFIED

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}


DSP_RTN_CODE DecBlkModel::decomposeCoupling(
		int size,                    /**< [in] size of subproblem subset */
		int * subprobs,              /**< [in] subset of subproblems */
		CoinPackedMatrix *& cpl_mat, /**< [out] coupling constraint matrix */
		int *& cpl_cols,             /**< [out] columns of cpl_mat involved in coupling rows */
		int & cpl_ncols              /**< [out] size of cpl_cols */) {

	BGN_TRY_CATCH

	/** retrieve the master matrix */
	const CoinPackedMatrix* mat = blk_->block(0)->getConstraintMatrix();

	/** coupling rows and columns */
	std::vector<int> coupling_rows;
	std::vector<int> coupling_cols;
	for (int s = 0; s < size; ++s) {
		DetBlock* sub = blk_->block(subprobs[s]+1);
		for (int i = 0; i < sub->getNumCouplingRows(); ++i)
			coupling_rows.push_back(sub->getCouplingRows()[i]);
		for (int i = 0; i < sub->getNumCouplingCols(); ++i)
			coupling_cols.push_back(sub->getCouplingCols()[i]);
	}

	/** TODO: These steps would better go into the CoinPackedVector class or a derived class. */
	/** erase duplicates */
	std::sort(coupling_rows.begin(), coupling_rows.end());
	coupling_rows.erase(
			std::unique(coupling_rows.begin(), coupling_rows.end()),
			coupling_rows.end());
	std::sort(coupling_cols.begin(), coupling_cols.end());
	coupling_cols.erase(
			std::unique(coupling_cols.begin(), coupling_cols.end()),
			coupling_cols.end());

	/** create a sub-matrix of the master matrix */
	cpl_mat->submatrixOf(*mat, coupling_rows.size(), &coupling_rows[0]);

	/** assign coupilng columns */
	cpl_ncols = coupling_cols.size();
	cpl_cols = new int [cpl_ncols];
	CoinCopyN(&coupling_cols[0], cpl_ncols, cpl_cols);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DecBlkModel::getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) {

	BGN_TRY_CATCH

	/** the master problem part */
	DetBlock* master = blk_->block(0);
	/** number of rows and columns */
	int nrows = master->getNumRows();
	int ncols = master->getNumCols();
	/** constraint matrix */
	mat = new CoinPackedMatrix(*(master->getConstraintMatrix()));
	/** make sure the matrix is row-wise */
	if (mat->isColOrdered()) mat->reverseOrdering();

	/** the subproblem part */
	int coloffset = 0;
	for (int i = 1; i < blk_->getNumBlocks(); ++i) {
		/** retrieve sub-block */
		DetBlock* sub = blk_->block(i);
		/** number of rows and columns */
		nrows += sub->getNumRows();
		ncols += sub->getNumCols();
		/** constraint matrix */
		const CoinPackedMatrix* submat = sub->getConstraintMatrix();
		const int* indices = submat->getIndices();
		const int* start = submat->getVectorStarts();
		const double* values = submat->getElements();
		/** add rows with column indices adjusted */
		for (int k = 0; k < submat->getNumRows(); ++k) {
			CoinPackedVector row;
			for (int j = 0; j < submat->getVectorSize(k); ++j) {
				if (indices[start[k]+j] < mat->getNumCols())
					row.insert(indices[start[k]+j], values[start[k]+j]);
				else
					row.insert(indices[start[k]+j] + coloffset, values[start[k]+j]);
			}
			mat->appendRow(row);
		}
		coloffset += sub->getNumCols() - mat->getNumCols();
		sub = NULL;
	}

	/** allocate memory for the outputs */
	clbd  = new double [ncols];
	cubd  = new double [ncols];
	ctype = new char [ncols];
	obj   = new double [ncols];
	rlbd  = new double [nrows];
	rubd  = new double [nrows];

	/** copy the others */
	CoinCopyN(master->getColLower(), master->getNumCols(), clbd);
	CoinCopyN(master->getColUpper(), master->getNumCols(), cubd);
	CoinCopyN(master->getCtype(), master->getNumCols(), ctype);
	CoinCopyN(master->getObj(), master->getNumCols(), obj);
	CoinCopyN(master->getRowLower(), master->getNumRows(), rlbd);
	CoinCopyN(master->getRowUpper(), master->getNumRows(), rubd);
	/** mark positions */
	int cpos = master->getNumCols();
	int rpos = master->getNumRows();
	/** for sub-blocks */
	for (int i = 1; i < blk_->getNumBlocks(); ++i) {
		DetBlock* sub = blk_->block(i);
		CoinCopyN(sub->getColLower(), sub->getNumCols(), clbd + cpos);
		CoinCopyN(sub->getColUpper(), sub->getNumCols(), cubd + cpos);
		CoinCopyN(sub->getCtype(), sub->getNumCols(), ctype + cpos);
		CoinCopyN(sub->getObj(), sub->getNumCols(), obj + cpos);
		CoinCopyN(sub->getRowLower(), sub->getNumRows(), rlbd + rpos);
		CoinCopyN(sub->getRowUpper(), sub->getNumRows(), rubd + rpos);
		cpos += sub->getNumCols();
		rpos += sub->getNumRows();
		sub = NULL;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE DecBlkModel::getFullModel(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		CoinPackedMatrix *& qobj,/**< [out] quadratic objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */) {

	BGN_TRY_CATCH

	/** the master problem part */
	DetBlock* master = blk_->block(0);
	/** number of rows and columns */
	int nrows = master->getNumRows();
	int ncols = master->getNumCols();
	/** constraint matrix */
	mat = new CoinPackedMatrix(*(master->getConstraintMatrix()));
	/** make sure the matrix is row-wise */
	if (mat->isColOrdered()) mat->reverseOrdering();

	/** master quadratic objective */
	/*
	CoinPackedMatrix *qmaster=NULL;
	if (master->getQuadraticObjectiveMatrix() != NULL){
		qmaster = new CoinPackedMatrix(*(master->getQuadraticObjectiveMatrix()));
		if (qmaster->isColOrdered()) qmaster->reverseOrdering();
		cout << "master qobj is not empty" <<endl;
	}
	*/

	/*		
	vector<int> qrowIndices;
	vector<int> qcolIndices;
	vector<double> qelements;
	*/

	/*
	int rownum=0;
	int pos=0;
	if (qmaster !=NULL){
		const CoinBigIndex * qstarts=qmaster->getVectorStarts();
		for (int i = 0; i < qmaster->getNumRows(); ++i)
		{
			int qlength = qmaster->getVectorSize(i);
			DSPdebugMessage("length of the %d th vector = %d", i, qlength);
			for (int j = 0; j < qlength; ++j) {
					qrowIndices.push_back(rownum);
					qcolIndices.push_back(qmaster->getIndices()[qstarts[i]+j]);
					if (rownum!=qmaster->getIndices()[qstarts[i]+j]){
						qelements.push_back(2*qmaster->getElements()[qstarts[i]+j]);
					}
					else{
						qelements.push_back(qmaster->getElements()[qstarts[i]+j]);
					}

					pos++;
			}	
			rownum++;
		}
	}*/
	
	/** the subproblem part */
	int coloffset = 0;
	for (int i = 1; i < blk_->getNumBlocks(); ++i) {
		/** retrieve sub-block */
		DetBlock* sub = blk_->block(i);
		/** number of rows and columns */
		nrows += sub->getNumRows();
		ncols += sub->getNumCols();
		/** constraint matrix */
		const CoinPackedMatrix* submat = sub->getConstraintMatrix();
		const int* indices = submat->getIndices();
		const int* start = submat->getVectorStarts();
		const double* values = submat->getElements();

		/** quadratic objective */
		/*
		const CoinPackedMatrix* subqobj = NULL;
		if (sub->getQuadraticObjectiveMatrix() != NULL){
			subqobj=new CoinPackedMatrix(*(sub->getQuadraticObjectiveMatrix()));
		}
		*/

		/** add rows with column indices adjusted */
		for (int k = 0; k < submat->getNumRows(); ++k) {
			CoinPackedVector row;
			for (int j = 0; j < submat->getVectorSize(k); ++j) {
				if (indices[start[k]+j] < mat->getNumCols())
					row.insert(indices[start[k]+j], values[start[k]+j]);
				else
					row.insert(indices[start[k]+j] + coloffset, values[start[k]+j]);
			}
			mat->appendRow(row);
		}

		/*
		if (subqobj != NULL){
			int rownum=0;
			for (int k = 0; k < subqobj->getNumRows(); ++k){
				CoinBigIndex qstart = subqobj->getVectorStarts()[k];
				for (int j=0; j < subqobj->getVectorSize(k); j++){
					if (subqobj->getIndices()[qstart+j]<mat->getNumCols()){
						qcolIndices.push_back(subqobj->getIndices()[qstart+j]);
					}
					else{
						qcolIndices.push_back(subqobj->getIndices()[qstart+j]+coloffset);
					}
					if(rownum<mat->getNumCols()){
						qrowIndices.push_back(rownum);
					}
					else{
						qrowIndices.push_back(rownum+coloffset);
					}
					if (rownum!=subqobj->getIndices()[qstart+j]){
						qelements.push_back(2 * (subqobj->getElements()[qstart + j]));
					}
					else{
						qelements.push_back((subqobj->getElements()[qstart + j]));
					}
				}
				rownum++;
			}
		}
		*/
		coloffset += sub->getNumCols() - mat->getNumCols();
		sub = NULL;
	}

	/** allocate memory for the outputs */
	clbd  = new double [ncols];
	cubd  = new double [ncols];
	ctype = new char [ncols];
	obj   = new double [ncols];
	rlbd  = new double [nrows];
	rubd  = new double [nrows];

	/** copy the others */
	CoinCopyN(master->getColLower(), master->getNumCols(), clbd);
	CoinCopyN(master->getColUpper(), master->getNumCols(), cubd);
	CoinCopyN(master->getCtype(), master->getNumCols(), ctype);
	CoinCopyN(master->getObj(), master->getNumCols(), obj);
	CoinCopyN(master->getRowLower(), master->getNumRows(), rlbd);
	CoinCopyN(master->getRowUpper(), master->getNumRows(), rubd);
	/** mark positions */
	int cpos = master->getNumCols();
	int rpos = master->getNumRows();
	/** for sub-blocks */
	for (int i = 1; i < blk_->getNumBlocks(); ++i) {
		DetBlock* sub = blk_->block(i);
		CoinCopyN(sub->getColLower(), sub->getNumCols(), clbd + cpos);
		CoinCopyN(sub->getColUpper(), sub->getNumCols(), cubd + cpos);
		CoinCopyN(sub->getCtype(), sub->getNumCols(), ctype + cpos);
		CoinCopyN(sub->getObj(), sub->getNumCols(), obj + cpos);
		CoinCopyN(sub->getRowLower(), sub->getNumRows(), rlbd + rpos);
		CoinCopyN(sub->getRowUpper(), sub->getNumRows(), rubd + rpos);
		cpos += sub->getNumCols();
		rpos += sub->getNumRows();
		sub = NULL;
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

void DecBlkModel::__printData() {
	for (int s = 0; s < blk_->getNumBlocks(); ++s) {
		printf("### Block Id: %d ###\n", s);
		blk_->block(s)->__printData();
	}
}
