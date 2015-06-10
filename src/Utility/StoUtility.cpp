/*
 * StoUtility.cpp
 *
 *  Created on: Feb 2, 2015
 *      Author: kibaekkim
 */

#include "math.h"
#include "assert.h"

/** DSP */
#include "Utility/StoUtility.h"
#include "Utility/StoMacros.h"
#include "Utility/StoMessage.h"

/** check whether solution is duplicate or not; return NULL if duplicate */
bool duplicateVector(
		CoinPackedVector * vec,
		vector<CoinPackedVector*> vecs)
{
	bool dup = false;

	/** number of saved solutions */
	int num = vecs.size();
	for (int i = num - 1; i >= 0; --i)
	{
		if (vec->getNumElements() != vecs[i]->getNumElements() ||
				vec->getMinIndex() != vecs[i]->getMinIndex() ||
			vec->getMaxIndex() != vecs[i]->getMaxIndex() ||
			fabs(vec->infNorm() - vecs[i]->infNorm()) > 1.0e-8 ||
			fabs(vec->oneNorm() - vecs[i]->oneNorm()) > 1.0e-8 ||
			fabs(vec->sum() - vecs[i]->sum()) > 1.0e-8 ||
			fabs(vec->twoNorm() - vecs[i]->twoNorm()) > 1.0e-8)
			continue;
		if (vec->isEquivalent(*vecs[i]))
		{
			dup = true;
			break;
		}
	}

	return dup;
}

void MPIscatterCoinPackedVectors(
		MPI::Intracomm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out)
{
	int comm_size = comm.Get_size();
	int   number_of_vectors  = vecs_in.size(); /**< number of vectors in process */
	int * numbers_of_vectors = NULL;           /**< number of vectors per process */
	int   total_number_of_vectors = 0;         /**< number of vectors in comm world */
	int   number_of_elements = 0;     /**< number of elements in process */
	int * numbers_of_elements = NULL; /**< number of elements per vector in comm world */
	int   length_of_vectors = 0;      /**< number of elements in comm world */
	int * indices = NULL;
	double * values = NULL;

	/** receive buffer specific */
	int * rcounts = new int [comm_size]; /**< number of elements that are to be received from each process */
	int * displs  = new int [comm_size]; /**< Entry i specifies the displacement at which to place the incoming data from peocess i */

	/** all gather number of vectors */
	{
		int sendbuf[1] = {number_of_vectors};
		numbers_of_vectors = new int [comm_size];
		/** communicate */
		comm.Allgather(sendbuf, 1, MPI::INT, numbers_of_vectors, 1, MPI::INT);
	}

	/** all gather numbers of elements */
	{
		int * sendbuf = new int [number_of_vectors];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			DSPdebugMessage("-> Process %d: vecs_in[%d]->getNumElements() %d\n",
					comm.Get_rank(), i, vecs_in[i]->getNumElements());
			sendbuf[i] = vecs_in[i]->getNumElements();
		}
		for (int i = 0; i < comm_size; ++i)
		{
			rcounts[i] = numbers_of_vectors[i];
			displs[i] = i == 0 ? 0 : displs[i-1] + rcounts[i-1];
			/** get total size of input vectors */
			total_number_of_vectors += rcounts[i];
		}
		numbers_of_elements = new int [total_number_of_vectors];
		/** communicate */
		comm.Allgatherv(sendbuf, number_of_vectors, MPI::INT,
				numbers_of_elements, rcounts, displs, MPI::INT);
		/** free send buffer */
		FREE_ARRAY_PTR(sendbuf);
	}

	/** get number of elements per process */
	for (int i = 0; i < number_of_vectors; ++i)
		number_of_elements += vecs_in[i]->getNumElements();

	/** all gather indices */
	{
		int * sendbuf = new int [number_of_elements];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinCopyN(vecs_in[i]->getIndices(), vecs_in[i]->getNumElements(), sendbuf);
			sendbuf += vecs_in[i]->getNumElements();
		}
		sendbuf -= number_of_elements;
		for (int i = 0, k = 0; i < comm_size; ++i)
		{
			rcounts[i] = 0;
			displs[i] = 0;
			for (int j = 0; j < numbers_of_vectors[i]; ++j)
				rcounts[i] += numbers_of_elements[k++];
			displs[i] = i == 0 ? 0 : displs[i-1] + rcounts[i-1];
			/** get length of input vectors */
			length_of_vectors += rcounts[i];
		}
		indices = new int [length_of_vectors];
		/** communicate */
		comm.Allgatherv(sendbuf, number_of_elements, MPI::INT,
				indices, rcounts, displs, MPI::INT);
		/** free send buffer */
		FREE_ARRAY_PTR(sendbuf);
	}

	/** all gather values */
	{
		double * sendbuf = new double [number_of_elements];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinCopyN(vecs_in[i]->getElements(), vecs_in[i]->getNumElements(), sendbuf);
			sendbuf += vecs_in[i]->getNumElements();
		}
		sendbuf -= number_of_elements;
		values = new double [length_of_vectors];
		/** communicate */
		comm.Allgatherv(sendbuf, number_of_elements, MPI::DOUBLE,
				values, rcounts, displs, MPI::DOUBLE);
		/** free send buffer */
		FREE_ARRAY_PTR(sendbuf);
	}

	/** re-construct output vectors */
	{
		DSPdebugMessage("Rank %d: total_number_of_vectors %d\n",
				comm.Get_rank(), total_number_of_vectors);
		for (int i = 0; i < total_number_of_vectors; ++i)
		{
			CoinPackedVector * vec = new CoinPackedVector(
					numbers_of_elements[i], indices, values);
			indices += numbers_of_elements[i];
			values  += numbers_of_elements[i];
			vecs_out.push_back(vec);
		}
		indices -= length_of_vectors;
		values  -= length_of_vectors;
	}

	/** free memory */
	FREE_ARRAY_PTR(numbers_of_vectors);
	FREE_ARRAY_PTR(numbers_of_elements);
	FREE_ARRAY_PTR(indices);
	FREE_ARRAY_PTR(values);
	FREE_ARRAY_PTR(rcounts);
	FREE_ARRAY_PTR(displs);
}

/** MPI scatter function for OsiCuts type. */
void MPIscatterOsiCuts(
		MPI::Intracomm comm,
		OsiCuts cuts_in,
		OsiCuts & cuts_out)
{
	/** input: cut elements */
	int ncuts = cuts_in.sizeCuts(); /**< number of cuts in this processor */
	vector<CoinPackedVector*> rows_in;
	double * lhs_in = NULL;
	double * rhs_in = NULL;

	/** output: cut elements */
	int ncuts_out = 0;
	vector<CoinPackedVector*> rows_out;
	double * lhs_out = NULL;
	double * rhs_out = NULL;

	/** MPI data */
	int comm_size = comm.Get_size(); /**< number of processors */
	int * ncutss = NULL;             /**< number of cuts for each processor */
	int * displs = NULL;

	/** allocate memory */
	lhs_in = new double [ncuts];
	rhs_in = new double [ncuts];
	ncutss = new int [ncuts];
	displs = new int [comm_size];

	/** parse cuts */
	for (int i = 0; i < ncuts; ++i)
	{
		OsiRowCut * rc = cuts_in.rowCutPtr(i);
		assert(rc);
		/** row vectors */
		rows_in.push_back(new CoinPackedVector(rc->row()));
		/** lhs */
		lhs_in[i] = rc->lb();
		/** rhs */
		rhs_in[i] = rc->ub();
	}

	/** synchronize row vectors */
	MPIscatterCoinPackedVectors(comm, rows_in, rows_out);

	/** all gather number of vectors */
	{
		int sendbuf[1] = {ncuts};
		ncutss = new int [comm_size];
		/** communicate */
		comm.Allgather(sendbuf, 1, MPI::INT, ncutss, 1, MPI::INT);
	}

	{
		for (int i = 0; i < comm_size; ++i)
		{
			displs[i] = i == 0 ? 0 : displs[i-1] + ncutss[i-1];
			/** get length of input vectors */
			ncuts_out += ncutss[i];
		}
		//printf(" --> Process %d: ncuts_out %d\n", comm.Get_rank(), ncuts_out);
		lhs_out = new double [ncuts_out];
		rhs_out = new double [ncuts_out];
		/** synchronize lhs */
		comm.Allgatherv(lhs_in, ncuts, MPI::DOUBLE,
				lhs_out, ncutss, displs, MPI::DOUBLE);
		/** synchronize rhs */
		comm.Allgatherv(rhs_in, ncuts, MPI::DOUBLE,
				rhs_out, ncutss, displs, MPI::DOUBLE);
	}

	/** recover cuts */
	for (int i = 0; i < ncuts_out; ++i)
	{
		OsiRowCut * rc = new OsiRowCut;
		rc->setRow(*rows_out[i]);
		rc->setLb(lhs_out[i]);
		rc->setUb(rhs_out[i]);
		cuts_out.insert(rc);
	}

	/** release memeory */
	for (unsigned int i = 0; i < rows_in.size(); ++i)
		FREE_PTR(rows_in[i]);
	FREE_ARRAY_PTR(lhs_in);
	FREE_ARRAY_PTR(rhs_in);
	for (unsigned int i = 0; i < rows_out.size(); ++i)
		FREE_PTR(rows_out[i]);
	FREE_ARRAY_PTR(lhs_out);
	FREE_ARRAY_PTR(rhs_out);
	FREE_ARRAY_PTR(ncutss);
	FREE_ARRAY_PTR(displs);
}

