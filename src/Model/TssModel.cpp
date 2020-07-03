/*
 * TssModel.cpp
 *
 *  Created on: Sep 23, 2014
 *      Author: kibaekkim
 */

#define DSP_DEBUG

#include "Utility/DspMessage.h"
#include "Model/TssModel.h"

TssModel::TssModel() :
StoModel()
{
	nstgs_ = 2;
}

/** copy constructor */
TssModel::TssModel(const TssModel & rhs) :
StoModel(rhs)
{
	nstgs_ = 2;;
}

TssModel::~TssModel()
{
	/** nothing to do */
}

/** set number of scenarios */
DSP_RTN_CODE TssModel::setNumberOfScenarios(int nscen)
{
	BGN_TRY_CATCH

	nscen_ = nscen;

	/** allocate memory */
	nrows_      = new int [nstgs_];
	ncols_      = new int [nstgs_];
	nints_      = new int [nstgs_];
	rstart_     = new int [nstgs_];
	cstart_     = new int [nstgs_];
	clbd_core_  = new double * [nstgs_];
	cubd_core_  = new double * [nstgs_];
	obj_core_   = new double * [nstgs_];
	qobj_core_  = new CoinPackedMatrix * [nstgs_];
	rlbd_core_  = new double * [nstgs_];
	rubd_core_  = new double * [nstgs_];
	ctype_core_ = new char * [nstgs_];
	prob_       = new double [nscen_];
	mat_scen_   = new CoinPackedMatrix * [nscen_];
	clbd_scen_  = new CoinPackedVector * [nscen_];
	cubd_scen_  = new CoinPackedVector * [nscen_];
	obj_scen_   = new CoinPackedVector * [nscen_];
	qobj_scen_  = new CoinPackedMatrix * [nscen_];
	rlbd_scen_  = new CoinPackedVector * [nscen_];
	rubd_scen_  = new CoinPackedVector * [nscen_];

	/** initialize memory */
	for (int s = 0; s < nstgs_; ++s)
	{
		clbd_core_[s]  = NULL;
		cubd_core_[s]  = NULL;
		ctype_core_[s] = NULL;
		obj_core_[s]   = NULL;
		qobj_core_[s]  = NULL;
		rlbd_core_[s]  = NULL;
		rubd_core_[s]  = NULL;
	}
	for (int s = 0; s < nscen_; ++s)
	{
		mat_scen_[s]  = NULL;
		clbd_scen_[s] = NULL;
		cubd_scen_[s] = NULL;
		obj_scen_[s]  = NULL;
		qobj_scen_[s] = NULL;
		rlbd_scen_[s] = NULL;
		rubd_scen_[s] = NULL;
	}
	CoinZeroN(prob_, nscen_);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** set dimensions */
DSP_RTN_CODE TssModel::setDimensions(
		const int ncols1, /**< number of first-stage columns */
		const int nrows1, /**< number of first-stage rows */
		const int ncols2, /**< number of second-stage columns */
		const int nrows2  /**< number of second-stage rows */)
{
	BGN_TRY_CATCH

	/** core problem dimension */
	nrows_core_ = nrows1 + nrows2;
	ncols_core_ = ncols1 + ncols2;
	DSPdebugMessage("nrows_core_ %d ncols_core_ %d\n", nrows_core_, ncols_core_);

	/** stage information */
	nrows_[0]  = nrows1;
	ncols_[0]  = ncols1;
	nints_[0]  = 0;
	rstart_[0] = 0;
	cstart_[0] = 0;
	nrows_[1]  = nrows2;
	ncols_[1]  = ncols2;
	nints_[1]  = 0;
	rstart_[1] = nrows1;
	cstart_[1] = ncols1;

	/** allocate memory */
	for (int s = 0; s < nstgs_; ++s)
	{
		clbd_core_[s]  = new double [ncols_[s]];
		cubd_core_[s]  = new double [ncols_[s]];
		obj_core_[s]   = new double [ncols_[s]];
		rlbd_core_[s]  = new double [nrows_[s]];
		rubd_core_[s]  = new double [nrows_[s]];
		ctype_core_[s] = new char [ncols_[s]];
	}
	rows_core_ = new CoinPackedVector * [nrows_core_];

	/** initialize memory */
	for (int i = 0; i < nrows_core_; ++i)
		rows_core_[i] = NULL;

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** load first-stage problem */
DSP_RTN_CODE TssModel::loadFirstStage(
		const CoinBigIndex * start, /**< start index for each row */
		const int *          index, /**< column indices */
		const double *       value, /**< constraint elements */
		const double *       clbd,  /**< column lower bounds */
		const double *       cubd,  /**< column upper bounds */
		const char *         ctype, /**< column types */
		const double *       obj,   /**< objective coefficients */
		const double *       rlbd,  /**< row lower bounds */
		const double *       rubd   /**< row upper bounds */)
{
	BGN_TRY_CATCH

	if (ncols_ == NULL || ncols_[0] <= 0)
	{
		printf("Error: invalid number of columns.\n");
		return DSP_RTN_ERR;
	}

	if (nrows_ == NULL || nrows_[0] < 0)
	{
		printf("Error: invalid number of rows.\n");
		return DSP_RTN_ERR;
	}

	/** allocate memory and values */
	CoinCopyN(clbd,  ncols_[0], clbd_core_[0]);
	CoinCopyN(cubd,  ncols_[0], cubd_core_[0]);
	CoinCopyN(ctype, ncols_[0], ctype_core_[0]);
	CoinCopyN(obj,   ncols_[0], obj_core_[0]);
	CoinCopyN(rlbd,  nrows_[0], rlbd_core_[0]);
	CoinCopyN(rubd,  nrows_[0], rubd_core_[0]);

	/** count number of integer variables */
	nints_core_ = 0;
	for (int j = 0; j < ncols_[0]; ++j)
	{
		if (ctype_core_[0][j] != 'C')
		{
			nints_[0]++;
			nints_core_++;
		}
	}

	/** construct core matrix rows */
	for (int i = 0; i < nrows_[0]; ++i)
	{
		rows_core_[i] = new CoinPackedVector(start[i+1] - start[i], index + start[i], value + start[i]);
	}

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE TssModel::loadFirstStage(
		const CoinBigIndex * start, /**< start index for each row */
		const int *          index, /**< column indices */
		const double *       value, /**< constraint elements */
		const double *       clbd,  /**< column lower bounds */
		const double *       cubd,  /**< column upper bounds */
		const char *         ctype, /**< column types */
		const double *       obj,   /**< objective coefficients */
		const int * 		 qobjrowindex, /**< quadratic objective row indices */
		const int *			 qobjcolindex, /**< quadratic objective column indices */
		const double *		 qobjvalue, /**< quadratic objective constraint elements value */
		const int 			 numq,  /**< number of quadratic terms */
		const double *       rlbd,  /**< row lower bounds */
		const double *       rubd   /**< row upper bounds */)
{
	BGN_TRY_CATCH

	/** load linear part */
	loadFirstStage(start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);


	/** load quadratic objective */
	qobj_core_[0] = new CoinPackedMatrix(false, qobjrowindex, qobjcolindex, qobjvalue, numq);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/** load second-stage problem */
DSP_RTN_CODE TssModel::loadSecondStage(
        const int s,               /**< scenario index */
        const double prob,         /**< probability */
        const CoinBigIndex *start, /**< start index for each row */
        const int *index,          /**< column indices */
        const double *value,       /**< constraint elements */
        const double *clbd,        /**< column lower bounds */
        const double *cubd,        /**< column upper bounds */
        const char *ctype,         /**< column types */
        const double *obj,         /**< objective coefficients */
        const double *rlbd,        /**< row lower bounds */
        const double *rubd         /**< row upper bounds */) {
#define FREE_MEMORY \
        FREE_ARRAY_PTR(len); \
        FREE_ARRAY_PTR(cind); \
        FREE_ARRAY_PTR(rind);

    int *len = NULL;
    int *cind = NULL;
    int *rind = NULL;

    BGN_TRY_CATCH

        if (ncols_ == NULL || ncols_[1] <= 0) {
            printf("Error: invalid number of columns for scenario %d.\n", s);
            return DSP_RTN_ERR;
        }

        if (nrows_ == NULL || nrows_[1] <= 0) {
            printf("Error: invalid number of rows for scenario %d.\n", s);
            return DSP_RTN_ERR;
        }

        if (nscen_ <= 0) {
            printf("Error: invalid number of scenarios.\n");
            return DSP_RTN_ERR;
        }

        if (s < 0 || s >= nscen_) {
            printf("Error: invalid scenario index.\n");
            return DSP_RTN_ERR;
        }

        /** allocate values */
        prob_[s] = prob;
        CoinCopyN(clbd, ncols_[1], clbd_core_[1]);
        CoinCopyN(cubd, ncols_[1], cubd_core_[1]);
        CoinCopyN(ctype, ncols_[1], ctype_core_[1]);
        CoinCopyN(obj, ncols_[1], obj_core_[1]);
        CoinCopyN(rlbd, nrows_[1], rlbd_core_[1]);
        CoinCopyN(rubd, nrows_[1], rubd_core_[1]);
        //PRINT_ARRAY_MSG(ncols_[1], ctype_core_[1], "ctype_core_[1]");

        /** count number of integer variables */
        if (s == 0) {
            for (int j = 0; j < ncols_[1]; ++j) {
                if (ctype_core_[1][j] != 'C') {
                    nints_[1]++;
                    nints_core_++;
                }
            }
        }

        /** construct core matrix rows */
        for (int i = 0; i < nrows_[1]; ++i) {
            if (rows_core_[rstart_[1] + i] == NULL) {
                //printf("creating rows_core_[%d]\n", rstart_[1]+i);
                rows_core_[rstart_[1] + i] = new CoinPackedVector(start[i + 1] - start[i], index + start[i], 0.0);
                //rows_core_[rstart_[1] + i] = new CoinPackedVector(start[i+1] - start[i], index + start[i], value + start[i]);
            } else {
                for (int j = start[i]; j < start[i + 1]; ++j) {
                    bool added = false;
                    if (rows_core_[rstart_[1] + i]->findIndex(index[j]) < 0) {
                        //printf(" added index %d element %e to rows_core_[%d]\n", index[j], value[j], rstart_[1]+i);
                        rows_core_[rstart_[1] + i]->insert(index[j], 0.);
                        added = true;
                    }
                    if (added) rows_core_[rstart_[1] + i]->sortIncrIndex();
                }
            }
        }

        /** add scenario mapping */
        scen2stg_.insert(std::pair<int, int>(s, 1));

        /** row length */
        len = new int[nrows_[1]];
        for (int i = 0; i < nrows_[1]; ++i)
            len[i] = start[i + 1] - start[i];
        cind = new int[ncols_[1]];
        rind = new int[nrows_[1]];
        CoinIotaN(cind, ncols_[1], cstart_[1]);
        CoinIotaN(rind, nrows_[1], rstart_[1]);

        /** allocate memory */
        mat_scen_[s] = new CoinPackedMatrix(false, ncols_[0] + ncols_[1], nrows_[1], start[nrows_[1]], value, index,
                                            start, len);
        clbd_scen_[s] = new CoinPackedVector(ncols_[1], cind, clbd);
        cubd_scen_[s] = new CoinPackedVector(ncols_[1], cind, cubd);
        obj_scen_[s] = new CoinPackedVector(ncols_[1], cind, obj);
        rlbd_scen_[s] = new CoinPackedVector(nrows_[1], rind, rlbd);
        rubd_scen_[s] = new CoinPackedVector(nrows_[1], rind, rubd);
        DSPdebug(DspMessage::printArray(start[nrows_[1]], value));
        DSPdebug(mat_scen_[s]->verifyMtx(4));

    END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

    FREE_MEMORY;

    return DSP_RTN_OK;
}

DSP_RTN_CODE TssModel::loadSecondStage(
        const int s,               /**< scenario index */
        const double prob,         /**< probability */
        const CoinBigIndex *start, /**< start index for each row */
        const int *index,          /**< column indices */
        const double *value,       /**< constraint elements */
        const double *clbd,        /**< column lower bounds */
        const double *cubd,        /**< column upper bounds */
        const char *ctype,         /**< column types */
        const double *obj,         /**< objective coefficients */
		const int * qobjrowindex,  /**< quadratic objective row indices */
		const int *	qobjcolindex,  /**< quadratic objective column indices */
		const double * qobjvalue,  /**< quadratic objective constraint elements value */
        const int numq,			   /**< number of quadratic terms */
		const double *rlbd,        /**< row lower bounds */
        const double *rubd         /**< row upper bounds */) {
	
	loadSecondStage(s, prob, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);

	BGN_TRY_CATCH
	
	//not sure
	qobj_core_[1] = new CoinPackedMatrix(false, qobjrowindex, qobjcolindex, qobjvalue, numq);
	qobj_scen_[s] = new CoinPackedMatrix(false, qobjrowindex, qobjcolindex, qobjvalue, numq);
	DSPdebug(qobj_scen_[s]->verifyMtx(4));

	END_TRY_CATCH_RTN(;, DSP_RTN_ERR)
	return DSP_RTN_OK;
}

/**
 * This constructs a deterministic equivalent form.
 */
DSP_RTN_CODE TssModel::copyDeterministicEquivalent(
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		CoinPackedMatrix *&qobj, /**< [out] quadratic objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */)
{
	int * scen = new int [nscen_];
	CoinIotaN(scen, nscen_, 0);

	DSP_RTN_CODE rtn = decompose(nscen_, scen, 0, NULL, NULL, NULL,
			mat, clbd, cubd, ctype, obj, qobj, rlbd, rubd);

	/** free array */
	FREE_ARRAY_PTR(scen);

	return rtn;
}

/**
 * Create a DetModel representing the deterministic equivalent of this model.
 */
DSP_RTN_CODE TssModel::copyDeterministicEquivalent(
		DetModel *& det /**< [out] deterministic equivalent model */)
{
	CoinPackedMatrix * mat = NULL;
	double * clbd = NULL;
	double * cubd = NULL;
	char   * ctype = NULL;
	double * obj = NULL;
	CoinPackedMatrix * qobj = NULL;
	double * rlbd = NULL;
	double * rubd = NULL;

	DSP_RTN_CODE rtn = copyDeterministicEquivalent(mat, clbd, cubd, ctype, obj, qobj, rlbd, rubd);
	det = new DetModel(mat, clbd, cubd, ctype, obj, qobj, rlbd, rubd);

	return rtn;
}

/**
 * This routine decomposes the model based on inputs. If size = 0, then
 * this results in a standard Benders decomposition structure. If size = 1, then
 * this results in a dual decomposition structure.
 */
DSP_RTN_CODE TssModel::decompose(
		int size,                /**< [in] size of scenario subset */
		int * scen,              /**< [in] subset of scenarios */
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
		double *& rubd           /**< [out] row upper bounds */)
{
	assert(size >= 0);
	assert(nrows_);
	assert(ncols_);
	assert(naux >= 0);
	assert(mat == NULL);

	int s, i, j;
	int nrows = nrows_[0] + size * nrows_[1];
	int ncols = ncols_[0] + size * ncols_[1] + naux;
	DSPdebugMessage("nrows %d ncols %d\n", nrows, ncols);

	BGN_TRY_CATCH

	for (s = 0; s < size; ++s)
		if (mat_scen_[scen[s]] == NULL)
			throw "Invalid model data (mat_scen_)";

	/** allocate memory */
	clbd = new double [ncols];
	cubd = new double [ncols];
	ctype = new char [ncols];
	obj = new double [ncols];
	rlbd = new double [nrows];
	rubd = new double [nrows];

#ifdef DSP_TIMING
	double stime = CoinCpuTime();
#endif

	vector<int> rowIndices;
	vector<int> colIndices;
	vector<double> elements;
	// int * rowIndices = NULL;
	// int * colIndices = NULL;
	// double * elements = NULL;
	double * denserow = NULL;

	int nzcnt = 0;
	/** # of nonzeros in the first-stage matrix A */
	for (i = 0; i < nrows_[0]; ++i) {
		for (j = 0; j < rows_core_[i]->getNumElements(); ++j) {
			if (fabs(rows_core_[i]->getElements()[j]) > 1.e-10) //why doing this?
				nzcnt++;
		}
	}
	/** # of nonzeros in the technology and recourse matrices */
	for (s = 0; s < size; ++s)
	{
		if (fromSMPS_)
		{
			for (i = nrows_[0]; i < nrows_core_; ++i) {
				nzcnt += rows_core_[i]->getNumElements();
			}
		}
		else
			nzcnt += mat_scen_[scen[s]]->getNumElements();
	}
	DSPdebugMessage("nzcnt %d\n", nzcnt);

	/** allocate memory */
	rowIndices.reserve(nzcnt);
	colIndices.reserve(nzcnt);
	elements.reserve(nzcnt);
	// rowIndices = new int [nzcnt];
	// colIndices = new int [nzcnt];
	// elements = new double [nzcnt];
	if (fromSMPS_)
		denserow = new double [ncols_core_]; //why define the dense row?

	/** construction */
	int pos = 0;
	int rownum = 0;
	/** for matrix A */
	for (i = 0; i < nrows_[0]; ++i)
	{
		for (j = 0; j < rows_core_[i]->getNumElements(); ++j) {
			if (fabs(rows_core_[i]->getElements()[j]) > 1.e-10) {
				rowIndices.push_back(rownum);
				colIndices.push_back(rows_core_[i]->getIndices()[j]);
				elements.push_back(rows_core_[i]->getElements()[j]);
				pos++;
			}
		}
		rownum++;
	}
	for (s = 0; s < size; ++s)
	{
		for (i = nrows_[0]; i < nrows_core_; ++i)
		{
			CoinBigIndex start = mat_scen_[scen[s]]->getVectorStarts()[i - nrows_[0]];
			int length = mat_scen_[scen[s]]->getVectorSize(i - nrows_[0]);

			if (fromSMPS_)
			{
				/** create dense row */
				CoinZeroN(denserow, ncols_core_);
				for (int j = 0; j < rows_core_[i]->getNumElements(); ++j)
					denserow[rows_core_[i]->getIndices()[j]] = rows_core_[i]->getElements()[j];
				for (int j = 0; j < length; ++j)
					denserow[mat_scen_[scen[s]]->getIndices()[start + j]] = mat_scen_[scen[s]]->getElements()[start + j];

				/** create sparse row */
				length = 0;
				for (int j = 0; j < ncols_core_; ++j)
				{
					if (fabs(denserow[j]) > 1.e-10)
					{
						rowIndices.push_back(rownum);
						colIndices.push_back(j);
						elements.push_back(denserow[j]);
						length++;
					}
				}
			}
			else /** from Julia */
			{
				for (j = 0; j < length; ++j) {
					rowIndices.push_back(rownum);
					colIndices.push_back(mat_scen_[scen[s]]->getIndices()[start + j]);
					elements.push_back(mat_scen_[scen[s]]->getElements()[start + j]);
				}
			}
			shiftVecIndices(length, &colIndices[0] + pos, s * ncols_[1], cstart_[1]);
			DSPdebugMessage("shift vector indices from %d by %d from %d\n", pos, s * ncols_[1], cstart_[1]);
			/*for (int k = 0, l = 0; l < length; ++l)
			{
				if (fabs(elements[pos+l]) < 1.0e-10) continue;
				if (k > 0 && k % 4 == 0) printf("\n");
				printf("  [%5d,%4d,%4d] %+e", i, rowIndices[pos+l], mat_scen_[scen[s]]->getIndices()[start+l], mat_scen_[scen[s]]->getElements()[start+l]);
				//printf("  [%5d,%4d,%4d] %+e", i, rowIndices[pos+l], colIndices[pos+l], elements[pos+l]);
				k++;
			}
			printf("\n");*/
			pos += length;
			rownum++;
		}
	}
	nzcnt = rowIndices.size();
	DSPdebugMessage("rownum %d nzcnt %d pos %d\n", rownum, nzcnt, pos);
	assert(nzcnt == pos);

	mat = new CoinPackedMatrix(false, &rowIndices[0], &colIndices[0], &elements[0], nzcnt);
	mat->setDimensions(nrows, ncols);
	DSPdebug(mat->verifyMtx(4));

	/** free memory */
	rowIndices.clear();
	colIndices.clear();
	elements.clear();
	FREE_ARRAY_PTR(denserow);

#ifdef DSP_TIMING
	printf("construct matrix %f seconds.\n", CoinCpuTime() - stime);
#endif

	/** all the other model data */
	copyCoreColLower(clbd, 0);
	copyCoreColUpper(cubd, 0);
	copyCoreColType(ctype, 0);
	copyCoreObjective(obj, 0);
	copyCoreRowLower(rlbd, 0);
	copyCoreRowUpper(rubd, 0);
	for (s = 0; s < size; ++s)
	{
		int sind = scen[s];

		/** column lower bounds */
		DSPdebugMessage("combining column lower bounds\n");
		copyCoreColLower(clbd + ncols_[0] + s * ncols_[1], 1);
		combineRandColLower(clbd + ncols_[0] + s * ncols_[1], 1, sind);

		/** column upper bounds */
		DSPdebugMessage("combining column upper bounds\n");
		copyCoreColUpper(cubd + ncols_[0] + s * ncols_[1], 1);
		combineRandColUpper(cubd + ncols_[0] + s * ncols_[1], 1, sind);

		/** column types */
		DSPdebugMessage("combining column types\n");
		copyCoreColType(ctype + ncols_[0] + s * ncols_[1], 1);

		/** objective coefficients */
		DSPdebugMessage("combining objective coefficients\n");
		copyCoreObjective(obj + ncols_[0] + s * ncols_[1], 1);
		combineRandObjective(obj + ncols_[0] + s * ncols_[1], 1, sind);

		/** row lower bounds */
		DSPdebugMessage("combining row lower bounds\n");
		copyCoreRowLower(rlbd + nrows_[0] + s * nrows_[1], 1);
		combineRandRowLower(rlbd + nrows_[0] + s * nrows_[1], 1, sind);

		/** row upper bounds */
		DSPdebugMessage("combining row upper bounds\n");
		copyCoreRowUpper(rubd + nrows_[0] + s * nrows_[1], 1);
		combineRandRowUpper(rubd + nrows_[0] + s * nrows_[1], 1, sind);
	}

	/** auxiliary columns */
	CoinCopyN(clbd_aux, naux, clbd + ncols - naux);
	CoinCopyN(cubd_aux, naux, cubd + ncols - naux);
	CoinCopyN(obj_aux, naux, obj + ncols - naux);
	CoinFillN(ctype + ncols - naux, naux, 'C');

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE TssModel::decompose(
		int size,                /**< [in] size of scenario subset */
		int * scen,              /**< [in] subset of scenarios */
		int naux,                /**< [in] number of auxiliary columns */
		double * clbd_aux,       /**< [in] lower bounds for auxiliary columns */
		double * cubd_aux,       /**< [in] upper bounds for auxiliary columns */
		double * obj_aux,        /**< [in] objective coefficients for auxiliary columns */
		//CoinPackedMatrix * qobj_aux, /**< [in] quadratic objective coefficients for auxiliary columns */
		CoinPackedMatrix *& mat, /**< [out] constraint matrix */
		double *& clbd,          /**< [out] column lower bounds */
		double *& cubd,          /**< [out] column upper bounds */
		char   *& ctype,         /**< [out] column types */
		double *& obj,           /**< [out] objective coefficients */
		CoinPackedMatrix * &qobj,/**< [out] quadratic objective coefficients */
		double *& rlbd,          /**< [out] row lower bounds */
		double *& rubd           /**< [out] row upper bounds */)
{
	decompose(size, scen, naux, clbd_aux, cubd_aux, obj_aux, mat, clbd, cubd, ctype, obj, rlbd, rubd);
	
	assert(qobj == NULL);
	int s;

	BGN_TRY_CATCH

	for (s = 0; s < size; ++s)
		if (qobj_scen_[scen[s]] == NULL)
			throw "Invalid model data (qobj_scen_)";
	
	/** # of quadratic objective elements in the first stage */
	int nzcnt = qobj_core_[0]->getNumElements();
	int nqrow = qobj_core_[0]->getMajorDim();

	/** # of quadratic terms in the second stage */
	/*
	for (s=0; s < size; s++){
		//TODO: from SMPS?

	}
	*/

	vector<int> rowIndices;
	vector<int> colIndices;
	vector<double> elements;
	int pos = 0;
	int rownum = 0;
	for (s=0; s<size; ++s){
		for (int i = 0; i<qobj_scen_[scen[s]]->getMajorDim() ;++i){
			CoinBigIndex start = qobj_scen_[scen[s]]->getVectorStarts()[i];
			int length = qobj_scen_[scen[s]]->getVectorSize(i);
			for (int j=0; j<length; ++j){
				rowIndices.push_back(rownum);
				colIndices.push_back(qobj_scen_[scen[s]]->getIndices()[start + j]);
				elements.push_back(qobj_scen_[scen[s]]->getElements()[start + j]);
			}
			shiftVecIndices(length, &colIndices[0] + pos, s * ncols_[1], cstart_[1]);
			DSPdebugMessage("shift vector in qobj indices from %d by %d from %d\n", pos, s * ncols_[1], cstart_[1]);
			pos += length;
			rownum++;
		}
	}
	nzcnt = rowIndices.size();
	DSPdebugMessage("qobj rownum %d nzcnt %d pos %d\n", rownum, nzcnt, pos);
	assert(nzcnt == pos);
	qobj = new CoinPackedMatrix(false, &rowIndices[0], &colIndices[0], &elements[0], nzcnt);
	//qobj->setDimensions(nrows, ncols);
	DSPdebug(qobj->verifyMtx(4));

	/** free memory */
	rowIndices.clear();
	colIndices.clear();
	elements.clear();

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

/**
 * This copies recourse problem structure for a given scenario index.
 */
DSP_RTN_CODE TssModel::copyRecoProb(
		int scen,                     /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech, /**< [out] technology matrix */
		CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix */
		double *& clbd_reco,          /**< [out] column lower bounds */
		double *& cubd_reco,          /**< [out] column upper bounds */
		char   *& ctype_reco,         /**< [out] column types */
		double *& obj_reco,           /**< [out] objective coefficients */
		double *& rlbd_reco,          /**< [out] row lower bounds */
		double *& rubd_reco           /**< [out] row upper bounds */)
{
	assert(scen >= 0 && scen < nscen_);
	assert(nstgs_ > 1);
	assert(nrows_ != NULL);
	assert(ncols_ != NULL);
	assert(nrows_[0] >= 0);
	assert(ncols_[0] >= 0);
	assert(nrows_[1] >= 0);
	assert(ncols_[1] >= 0);

	int i;
	CoinPackedVector * row = NULL;

	BGN_TRY_CATCH

	/** allocate memory */
	mat_tech   = new CoinPackedMatrix(false, 0, 0);
	mat_reco   = new CoinPackedMatrix(false, 0, 0);
	clbd_reco  = new double [ncols_[1]];
	cubd_reco  = new double [ncols_[1]];
	ctype_reco = new char [ncols_[1]];
	rlbd_reco  = new double [nrows_[1]];
	rubd_reco  = new double [nrows_[1]];

	/** matrices */
	mat_tech->setDimensions(0, ncols_[0]);
	mat_reco->setDimensions(0, ncols_[1]);
	for (i = rstart_[1]; i < nrows_core_; ++i)
	{
		/** add row to tech */
		row = this->splitCoreRowVec(i, 0);
		combineRandRowVec(row, i - rstart_[1], 0, scen);
		mat_tech->appendRow(*row);
//		PRINT_COIN_PACKED_VECTOR_MSG((*row),"row_tech")
		FREE_PTR(row)

		/** add row to reco */
		row = this->splitCoreRowVec(i, 1);
		combineRandRowVec(row, i - rstart_[1], 1, scen);
		mat_reco->appendRow(*row);
		//PRINT_COIN_PACKED_VECTOR_MSG((*row),"row_reco")
		FREE_PTR(row)
	}

	/** column lower bounds */
	copyCoreColLower(clbd_reco, 1);
	combineRandColLower(clbd_reco, 1, scen);

	/** column upper bounds */
	copyCoreColUpper(cubd_reco, 1);
	combineRandColUpper(cubd_reco, 1, scen);

	/** column types */
	copyCoreColType(ctype_reco, 1);

	/** objective coefficients */
	copyRecoObj(scen, obj_reco, true);

	/** row lower bounds */
	copyCoreRowLower(rlbd_reco, 1);
	combineRandRowLower(rlbd_reco, 1, scen);

	/** row upper bounds */
	copyCoreRowUpper(rubd_reco, 1);
	combineRandRowUpper(rubd_reco, 1, scen);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE TssModel::copyRecoProb(
		int scen,                     /**< [in] scenario index */
		CoinPackedMatrix *& mat_tech, /**< [out] technology matrix */
		CoinPackedMatrix *& mat_reco, /**< [out] recourse matrix */
		double *& clbd_reco,          /**< [out] column lower bounds */
		double *& cubd_reco,          /**< [out] column upper bounds */
		char   *& ctype_reco,         /**< [out] column types */
		double *& obj_reco,           /**< [out] objective coefficients */
		CoinPackedMatrix *&qobj_reco, /**< [out] quadratic objective coefficients */
		double *& rlbd_reco,          /**< [out] row lower bounds */
		double *& rubd_reco           /**< [out] row upper bounds */)
{
	copyRecoProb(scen, mat_tech, mat_reco, clbd_reco, cubd_reco, ctype_reco, obj_reco, rlbd_reco, rubd_reco);

	BGN_TRY_CATCH

	copyRecoQuadraticObj(scen, qobj_reco, true);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE TssModel::copyRecoObj(int scen, double *& obj_reco, bool adjustProbability) {

	assert(ncols_[1] >= 0);

	BGN_TRY_CATCH

	/** allocate memory */
	obj_reco = new double [ncols_[1]];

	/** objective coefficients */
	copyCoreObjective(obj_reco, 1);
	combineRandObjective(obj_reco, 1, scen, adjustProbability);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;

}

DSP_RTN_CODE TssModel::copyRecoQuadraticObj(int scen, CoinPackedMatrix *& qobj_reco, bool adjustProbability) {

	BGN_TRY_CATCH

	copyCoreQuadrativeObjective(qobj_reco, 1);
	combineRandQuadraticObjective(qobj_reco, 1, scen, adjustProbability);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;

}