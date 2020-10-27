/*
 * TssModel.cpp
 *
 *  Created on: Sep 23, 2014
 *      Author: kibaekkim
 */

// #define DSP_DEBUG

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
	qobj_core_[0] = NULL;

	/** count number of integer variables */
	nints_core_ = 0;
	for (int j = 0; j < ncols_[0]; ++j)
	{
		if (ctype_core_[0][j] != 'C')
		{
			nints_[0]++;
			nints_core_++;

			/** set bounds for binary variables */
			if (ctype_core_[0][j] == 'B') {
				clbd_core_[0][j] = 0.0;
				cubd_core_[0][j] = 1.0;
			}
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
	if (qobjrowindex == NULL||qobjcolindex == NULL||qobjvalue == NULL || numq == 0){
		loadFirstStage(start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);
		return DSP_RTN_OK;
	}
	else{
		isQCQP_=1;
		loadFirstStage(start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);
		/** load quadratic objective */
		qobj_core_[0] = new CoinPackedMatrix(false, qobjrowindex, qobjcolindex, qobjvalue, numq);
		qobj_core_[0]->setDimensions(ncols_[0], ncols_[0]);
		DSPdebugMessage("get number of elements in qobjcore_[0] = %d\n", qobj_core_[0]->getNumElements());
		//PRINT_ARRAY_MSG(3, qobj_core_[0]->getVector(0), "the first row of qobj_core_[0]");
		//PRINT_ARRAY_MSG(3, qobj_core_[0]->getVector(1), "the second row of qobj_core_[0]");
	}
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

					/** set bounds for binary variables */
					if (ctype_core_[1][j] == 'B') {
						clbd_core_[1][j] = 0.0;
						cubd_core_[1][j] = 1.0;
					}
                }
            }
        }
		DSPdebugMessage("ncols_[1] = %d\n", ncols_[1]);
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
                        printf(" added index %d element %e to rows_core_[%d]\n", index[j], value[j], rstart_[1]+i);
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
        clbd_scen_[s] = new CoinPackedVector(ncols_[1], cind, clbd_core_[1]);
        cubd_scen_[s] = new CoinPackedVector(ncols_[1], cind, cubd_core_[1]);
        obj_scen_[s] = new CoinPackedVector(ncols_[1], cind, obj_core_[1]);
        rlbd_scen_[s] = new CoinPackedVector(nrows_[1], rind, rlbd_core_[1]);
        rubd_scen_[s] = new CoinPackedVector(nrows_[1], rind, rubd_core_[1]);
        DSPdebug(DspMessage::printArray(start[nrows_[1]], value));
        DSPdebug(mat_scen_[s]->verifyMtx(4));

    END_TRY_CATCH_RTN(FREE_MEMORY, DSP_RTN_ERR)

    FREE_MEMORY;

    return DSP_RTN_OK;
}

DSP_RTN_CODE TssModel::loadSecondStage(
			const int            s,     /**< scenario index */
			const double         prob,  /**< probability */
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
			const CoinBigIndex 	 qnum,  /**< number of quadratic terms */
			const double *       rlbd,  /**< row lower bounds */
			const double *       rubd   /**< row upper bounds */)
{
	BGN_TRY_CATCH
	if (qobjrowindex == NULL||qobjcolindex == NULL||qobjvalue == NULL || qnum == 0){
		loadSecondStage(s, prob, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);
	}
	else{
		isQCQP_=1;
		loadSecondStage(s, prob, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd);

		/** allocate memory for qobj_core_[1] */
		qobj_core_[1] = new CoinPackedMatrix(false, qobjrowindex, qobjcolindex, qobjvalue, qnum);
		qobj_core_[1]->setDimensions(ncols_[0] + ncols_[1], ncols_[0] + ncols_[1]);

		/** allocate memory for qobj_scen_ */
		qobj_scen_[s] = new CoinPackedMatrix(false, qobjrowindex, qobjcolindex, qobjvalue, qnum);
		qobj_scen_[s]->setDimensions(ncols_[0] + ncols_[1], ncols_[0] + ncols_[1]);
		DSPdebug(qobj_scen_[s]->verifyMtx(4));
	}
	DSPdebugMessage("load second stage with quadratic objective\n");
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
};

DSP_RTN_CODE TssModel::copyRecoObj(int scen, double *& obj_reco, bool adjustProbability) {

	assert(ncols_[1] >= 0);

	BGN_TRY_CATCH

	/** allocate memory */
	obj_reco = new double [ncols_[1]];

	/** objective coefficients */
	copyCoreObjective(obj_reco, 1);
	//PRINT_ARRAY_MSG(ncols_[1],obj_reco,"obj_reco");
	combineRandObjective(obj_reco, 1, scen, adjustProbability);
	//PRINT_ARRAY_MSG(ncols_[1],obj_reco,"obj_reco");
	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)
 
	return DSP_RTN_OK;

}

DSP_RTN_CODE TssModel::copyRecoObj(int scen, double *& obj_reco, CoinPackedMatrix *& qobj_reco_coupling, CoinPackedMatrix *& qobj_reco_ncoupling, bool adjustProbability) {

	BGN_TRY_CATCH

	copyRecoObj(scen, obj_reco, adjustProbability);

	qobj_reco_coupling=new CoinPackedMatrix(false, 0, 0);
	qobj_reco_ncoupling=new CoinPackedMatrix(false, 0, 0);

	/** copy quadratic objective coefficience to qobj_reco */
	copyCoreQuadrativeObjective(qobj_reco_coupling, qobj_reco_ncoupling, 1);
	//PRINT_ARRAY_MSG(qobj_reco_ncoupling->getNumElements(), qobj_reco_ncoupling->getElements(), "elements in qobj_reco_ncoupling11");
	combineRandQuadraticObjective(qobj_reco_coupling, qobj_reco_ncoupling, 1, scen, adjustProbability);
	//PRINT_ARRAY_MSG(qobj_reco_ncoupling->getNumElements(), qobj_reco_ncoupling->getElements(), "elements in qobj_reco_ncoupling");

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;

}