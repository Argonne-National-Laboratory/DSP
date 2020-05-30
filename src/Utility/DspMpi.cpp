/*
 * DspUtility.cpp
 *
 *  Refactored on: April 8, 2016
 *     Created on: Feb 2, 2015
 *         Author: kibaekkim
 */

//#define DSP_DEBUG

#include "math.h"
#include "assert.h"

#include "Utility/DspMacros.h"
#include "Utility/DspMessage.h"
#include "Utility/DspMpi.h"

/** DSP */

/** get round-and-robin distribution of indices */
DSP_RTN_CODE distIndices(
		int num_indices,       /**< [in] number of indices */
		int comm_size,         /**< [in] number of processors */
		int comm_rank,         /**< [in] processor id */
		int comm_rank_start,   /**< [in] smallest processor id */
		vector<int> & assigned /**< [out] assigned indices */)
{
	BGN_TRY_CATCH

	assigned.clear();
	for (int s = comm_rank - comm_rank_start; s < num_indices; s += comm_size)
		assigned.push_back(s);

	END_TRY_CATCH_RTN(;,DSP_RTN_ERR)

	return DSP_RTN_OK;
}

DSP_RTN_CODE MPIgatherCoinPackedVectors(
		MPI_Comm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out)
{
	int   number_of_vectors  = vecs_in.size(); /**< number of vectors in process */
	int * numbers_of_vectors = NULL;           /**< number of vectors per process */
	int   total_number_of_vectors = 0;         /**< number of vectors in comm world */
	int   number_of_elements = 0;     /**< number of elements in process */
	int * numbers_of_elements = NULL; /**< number of elements per vector in comm world */
	int   length_of_vectors = 0;      /**< number of elements in comm world */
	int * indices = NULL;
	double * values = NULL;

	/** receive buffer specific */
	int * rcounts = NULL; /**< number of elements that are to be received from each process */
	int * displs  = NULL; /**< Entry i specifies the displacement at which to place the incoming data from peocess i */

	int comm_size, comm_rank;
	MPI_Comm_size(comm, &comm_size);
	MPI_Comm_rank(comm, &comm_rank);

	if (comm_rank == 0)
	{
		rcounts = new int [comm_size];
		displs  = new int [comm_size];
	}

	/** all gather number of vectors */
	{
		int sendbuf[1] = {number_of_vectors};
		numbers_of_vectors = new int [comm_size];
		/** communicate */
		MPI_Gather(sendbuf, 1, MPI_INT, numbers_of_vectors, 1, MPI_INT, 0, comm);
	}

	/** all gather numbers of elements */
	{
		int * sendbuf = new int [number_of_vectors];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			DSPdebugMessage("-> Process %d: vecs_in[%d]->getNumElements() %d\n",
					comm_rank, i, vecs_in[i]->getNumElements());
			sendbuf[i] = vecs_in[i]->getNumElements();
		}
		if (comm_rank == 0)
		{
			for (int i = 0; i < comm_size; ++i)
			{
				rcounts[i] = numbers_of_vectors[i];
				displs[i] = i == 0 ? 0 : displs[i-1] + rcounts[i-1];
				/** get total size of input vectors */
				total_number_of_vectors += rcounts[i];
			}
			numbers_of_elements = new int [total_number_of_vectors];
		}
		/** communicate */
		MPI_Gatherv(sendbuf, number_of_vectors, MPI_INT,
				numbers_of_elements, rcounts, displs, MPI_INT, 0, comm);
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
		if (comm_rank == 0)
		{
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
		}
		/** communicate */
		MPI_Gatherv(sendbuf, number_of_elements, MPI_INT,
				indices, rcounts, displs, MPI_INT, 0, comm);
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
		if (comm_rank == 0)
			values = new double [length_of_vectors];
		/** communicate */
		MPI_Gatherv(sendbuf, number_of_elements, MPI_DOUBLE,
				values, rcounts, displs, MPI_DOUBLE, 0, comm);
		/** free send buffer */
		FREE_ARRAY_PTR(sendbuf);
	}

	/** re-construct output vectors */
	if (comm_rank == 0)
	{
		DSPdebugMessage("Rank %d: total_number_of_vectors %d\n",
				comm_rank, total_number_of_vectors);
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

	return DSP_RTN_OK;
}

DSP_RTN_CODE MPIAllgatherCoinPackedVectors(
		MPI_Comm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out)
{
	int   number_of_vectors  = vecs_in.size(); /**< number of vectors in process */
	int * numbers_of_vectors = NULL;           /**< number of vectors per process */
	int   total_number_of_vectors = 0;         /**< number of vectors in comm world */
	int   number_of_elements = 0;     /**< number of elements in process */
	int * numbers_of_elements = NULL; /**< number of elements per vector in comm world */
	int   length_of_vectors = 0;      /**< number of elements in comm world */
	int * indices = NULL;
	double * values = NULL;

	int comm_rank, comm_size;
	MPI_Comm_size(comm, &comm_size);
	MPI_Comm_rank(comm, &comm_rank);

	/** receive buffer specific */
	int * rcounts = new int [comm_size]; /**< number of elements that are to be received from each process */
	int * displs  = new int [comm_size]; /**< Entry i specifies the displacement at which to place the incoming data from peocess i */

	/** all gather number of vectors */
	{
		int sendbuf[1] = {number_of_vectors};
		numbers_of_vectors = new int [comm_size];
		/** communicate */
		MPI_Allgather(sendbuf, 1, MPI_INT, numbers_of_vectors, 1, MPI_INT, comm);
	}

	/** all gather numbers of elements */
	{
		int * sendbuf = new int [number_of_vectors];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			DSPdebugMessage("-> Process %d: vecs_in[%d]->getNumElements() %d\n",
					comm_rank, i, vecs_in[i]->getNumElements());
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
		MPI_Allgatherv(sendbuf, number_of_vectors, MPI_INT,
				numbers_of_elements, rcounts, displs, MPI_INT, comm);
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
		MPI_Allgatherv(sendbuf, number_of_elements, MPI_INT,
				indices, rcounts, displs, MPI_INT, comm);
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
		MPI_Allgatherv(sendbuf, number_of_elements, MPI_DOUBLE,
				values, rcounts, displs, MPI_DOUBLE, comm);
		/** free send buffer */
		FREE_ARRAY_PTR(sendbuf);
	}

	/** re-construct output vectors */
	{
		DSPdebugMessage("Rank %d: total_number_of_vectors %d\n",
				comm_rank, total_number_of_vectors);
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

	return DSP_RTN_OK;
}

/** MPI gather function for OsiCuts type. */
DSP_RTN_CODE MPIgatherOsiCuts(
		MPI_Comm comm,
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
	int comm_size;       /**< number of processors */
	int comm_rank;       /**< rank of process */
	int * ncutss = NULL; /**< number of cuts for each processor */
	int * displs = NULL;

	MPI_Comm_size(comm, &comm_size);
	MPI_Comm_rank(comm, &comm_rank);

	/** allocate memory */
	lhs_in = new double [ncuts];
	rhs_in = new double [ncuts];
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

	/** gather row vectors */
	MPIgatherCoinPackedVectors(comm, rows_in, rows_out);

	/** all gather number of vectors */
	{
		if (comm_rank == 0)
			ncutss = new int [comm_size];
		/** communicate */
		MPI_Gather(&ncuts, 1, MPI_INT, ncutss, 1, MPI_INT, 0, comm);
	}

	{
		if (comm_rank == 0)
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
		}
		/** synchronize lhs */
		MPI_Gatherv(lhs_in, ncuts, MPI_DOUBLE,
				lhs_out, ncutss, displs, MPI_DOUBLE, 0, comm);
		/** synchronize rhs */
		MPI_Gatherv(rhs_in, ncuts, MPI_DOUBLE,
				rhs_out, ncutss, displs, MPI_DOUBLE, 0, comm);
	}

	/** recover cuts */
	if (comm_rank == 0)
	{
		for (int i = 0; i < ncuts_out; ++i)
		{
			OsiRowCut * rc = new OsiRowCut;
			rc->setRow(*rows_out[i]);
			rc->setLb(lhs_out[i]);
			rc->setUb(rhs_out[i]);
			cuts_out.insert(rc);
		}
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

	return DSP_RTN_OK;
}

/** MPI Allgather function for OsiCuts type. */
DSP_RTN_CODE MPIAllgatherOsiCuts(
		MPI_Comm comm,
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
	int comm_size;       /**< number of processors */
	int * ncutss = NULL; /**< number of cuts for each processor */
	int * displs = NULL;
	MPI_Comm_size(comm, &comm_size);

	/** allocate memory */
	lhs_in = new double [ncuts];
	rhs_in = new double [ncuts];
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
	MPIAllgatherCoinPackedVectors(comm, rows_in, rows_out);

	/** all gather number of vectors */
	{
		ncutss = new int [comm_size];
		/** communicate */
		MPI_Allgather(&ncuts, 1, MPI_INT, ncutss, 1, MPI_INT, comm);
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
		MPI_Allgatherv(lhs_in, ncuts, MPI_DOUBLE,
				lhs_out, ncutss, displs, MPI_DOUBLE, comm);
		/** synchronize rhs */
		MPI_Allgatherv(rhs_in, ncuts, MPI_DOUBLE,
				rhs_out, ncutss, displs, MPI_DOUBLE, comm);
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

	return DSP_RTN_OK;
}

/** MPI_Send for vector<CoinPackedVectors*> */
DSP_RTN_CODE MPIsendCoinPackedVectors(
		MPI_Comm comm,
		int to_rank,
		vector<CoinPackedVector*> vecs,
		int tag)
{
#define FREE_MEMORY          \
	FREE_ARRAY_PTR(sendibuf) \
	FREE_ARRAY_PTR(senddbuf)

	/** send buffers */
	int *    sendibuf = NULL; /**< integer type */
	double * senddbuf = NULL; /**< double type */

	BGN_TRY_CATCH

	int number_of_vectors  = vecs.size(); /**< number of vectors */
	int number_of_elements = 0;           /**< number of elements */

	int comm_rank;
	MPI_Comm_rank(comm, &comm_rank);

	/** send the number of vectors */
	MPI_Send(&number_of_vectors, 1, MPI_INT, to_rank, tag, comm);

	/** send the numbers of elements */
	sendibuf = new int [number_of_vectors];
	for (int i = 0; i < number_of_vectors; ++i)
	{
		DSPdebugMessage("-> Process %d: vecs_in[%d]->getNumElements() %d\n",
				comm_rank, i, vecs[i]->getNumElements());
		sendibuf[i] = vecs[i]->getNumElements();
	}
	MPI_Send(sendibuf, number_of_vectors, MPI_INT, to_rank, tag, comm);
	for (int i = 0; i < number_of_vectors; ++i)
		number_of_elements += sendibuf[i];
	FREE_ARRAY_PTR(sendibuf);

	/** send indices */
	sendibuf = new int [number_of_elements];
	for (int i = 0; i < number_of_vectors; ++i)
	{
		CoinCopyN(vecs[i]->getIndices(), vecs[i]->getNumElements(), sendibuf);
		sendibuf += vecs[i]->getNumElements();
	}
	sendibuf -= number_of_elements;
	MPI_Send(sendibuf, number_of_elements, MPI_INT, to_rank, tag, comm);
	FREE_ARRAY_PTR(sendibuf);

	/** send values */
	senddbuf = new double [number_of_elements];
	for (int i = 0; i < number_of_vectors; ++i)
	{
		CoinCopyN(vecs[i]->getElements(), vecs[i]->getNumElements(), senddbuf);
		senddbuf += vecs[i]->getNumElements();
	}
	senddbuf -= number_of_elements;
	MPI_Send(senddbuf, number_of_elements, MPI_DOUBLE, to_rank, tag, comm);
	FREE_ARRAY_PTR(senddbuf);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** MPI_Recv for vector<CoinPackedVectors*> */
DSP_RTN_CODE MPIrecvCoinPackedVectors(
		MPI_Comm comm,
		int from_rank,
		vector<CoinPackedVector*> & vecs,
		int tag)
{
#define FREE_MEMORY          \
	FREE_ARRAY_PTR(indices) \
	FREE_ARRAY_PTR(elements) \
	FREE_ARRAY_PTR(numbers_of_elements)

	MPI_Status status;

	/** receive buffers */
	int *    indices  = NULL;
	double * elements = NULL;
	int * numbers_of_elements = NULL;

	BGN_TRY_CATCH

	int number_of_vectors;  /**< number of vectors */
	int number_of_elements = 0; /**< number of elements */

	/** receive the number of vectors */
	MPI_Recv(&number_of_vectors, 1, MPI_INT, from_rank, tag, comm, &status);

	/** receive the numbers of elements */
	numbers_of_elements = new int [number_of_vectors];
	MPI_Recv(numbers_of_elements, number_of_vectors, MPI_INT, from_rank, tag, comm, &status);
	for (int i = 0; i < number_of_vectors; ++i)
		number_of_elements += numbers_of_elements[i];

	/** receive indices */
	indices = new int [number_of_elements];
	MPI_Recv(indices, number_of_elements, MPI_INT, from_rank, tag, comm, &status);

	/** receive values */
	elements = new double [number_of_elements];
	MPI_Recv(elements, number_of_elements, MPI_DOUBLE, from_rank, tag, comm, &status);

	/** construct vectors */
	for (int i = 0; i < number_of_vectors; ++i)
	{
		CoinPackedVector * vec = new CoinPackedVector(
				numbers_of_elements[i], indices, elements);
		indices += numbers_of_elements[i];
		elements  += numbers_of_elements[i];
		vecs.push_back(vec);
	}
	indices -= number_of_elements;
	elements  -= number_of_elements;

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** MPI_Send for OsiCuts */
DSP_RTN_CODE MPIsendOsiCuts(
		MPI_Comm comm,
		int to_rank,
		OsiCuts cuts,
		int tag)
{
#define FREE_MEMORY     \
	FREE_ARRAY_PTR(lhs) \
	FREE_ARRAY_PTR(rhs)

	double * lhs = NULL;
	double * rhs = NULL;

	BGN_TRY_CATCH

	/** input: cut elements */
	int ncuts = cuts.sizeCuts(); /**< number of cuts in this processor */
	vector<CoinPackedVector*> rows;

	/** allocate memory */
	lhs = new double [ncuts];
	rhs = new double [ncuts];

	/** parse cuts */
	for (int i = 0; i < ncuts; ++i)
	{
		OsiRowCut * rc = cuts.rowCutPtr(i);
		assert(rc);
		/** row vectors */
		rows.push_back(new CoinPackedVector(rc->row()));
		/** lhs */
		lhs[i] = rc->lb();
		/** rhs */
		rhs[i] = rc->ub();
	}

	/** send row vectors */
	MPIsendCoinPackedVectors(comm, to_rank, rows, tag);

	/** send lhs */
	MPI_Send(lhs, ncuts, MPI_DOUBLE, to_rank, tag, comm);

	/** send rhs */
	MPI_Send(rhs, ncuts, MPI_DOUBLE, to_rank, tag, comm);

	/** free vector */
	for (unsigned i = 0; i < rows.size(); ++i)
		FREE_PTR(rows[i]);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** MPI_Recv for OsiCuts */
DSP_RTN_CODE MPIrecvOsiCuts(
		MPI_Comm comm,
		int from_rank,
		OsiCuts & cuts,
		int tag)
{
#define FREE_MEMORY     \
	FREE_ARRAY_PTR(lhs) \
	FREE_ARRAY_PTR(rhs)

	MPI_Status status;
	int ncuts;
	vector<CoinPackedVector*> rows;
	double * lhs = NULL;
	double * rhs = NULL;

	BGN_TRY_CATCH

	/** receive row vectors */
	MPIrecvCoinPackedVectors(comm, from_rank, rows, tag);

	ncuts = rows.size();
	lhs = new double [ncuts];
	rhs = new double [ncuts];

	/** receive lhs */
	MPI_Recv(lhs, ncuts, MPI_DOUBLE, from_rank, tag, comm, &status);

	/** receive rhs */
	MPI_Recv(rhs, ncuts, MPI_DOUBLE, from_rank, tag, comm, &status);

	/** construct cuts */
	for (int i = 0; i < ncuts; ++i)
	{
		OsiRowCut rc;
		rc.setRow(rows[i]);
		rc.setLb(lhs[i]);
		rc.setUb(rhs[i]);
		cuts.insert(rc);
	}

	/** free vector */
	for (unsigned i = 0; i < rows.size(); ++i)
		FREE_PTR(rows[i]);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** MPI_Scatter function for vector<CoinPackedVectors*> */
DSP_RTN_CODE MPIscatterCoinPackedVectors(
		MPI_Comm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(numbers_of_vectors)  \
	FREE_ARRAY_PTR(numbers_of_elements) \
	FREE_ARRAY_PTR(indices)             \
	FREE_ARRAY_PTR(values)

	int   number_of_vectors;          /**< number of vectors in process */
	int * numbers_of_vectors = NULL;  /**< number of vectors per process */
	int   number_of_elements;         /**< number of elements in process */
	int * numbers_of_elements = NULL; /**< number of elements per vector in comm world */
	int * indices = NULL;
	double * values = NULL;

	BGN_TRY_CATCH

	int comm_rank;
	MPI_Comm_rank(comm, &comm_rank);

	/** scatter the number of vectors */
	if (comm_rank == 0)
	{
		number_of_vectors = vecs_in.size();
		MPI_Scatter(&number_of_vectors, 1, MPI_INT, NULL, 0, MPI_INT, 0, comm);
	}
	else
		MPI_Scatter(NULL, 0, MPI_INT, &number_of_vectors, 1, MPI_INT, 0, comm);

	/** scatter the numbers of elements */
	if (comm_rank == 0)
	{
		numbers_of_elements = new int [number_of_vectors];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			DSPdebugMessage("-> Process %d: vecs_in[%d]->getNumElements() %d\n",
					comm_rank, i, vecs_in[i]->getNumElements());
			numbers_of_elements[i] = vecs_in[i]->getNumElements();
		}
		MPI_Scatter(numbers_of_elements, number_of_vectors, MPI_INT, NULL, 0, MPI_INT, 0, comm);
	}
	else
	{
		numbers_of_elements = new int [number_of_vectors];
		MPI_Scatter(NULL, 0, MPI_INT, numbers_of_elements, number_of_vectors, MPI_INT, 0, comm);
	}

	/** get number of elements per process */
	number_of_elements = 0;
	for (int i = 0; i < number_of_vectors; ++i)
		number_of_elements += numbers_of_elements[i];

	/** scatter indices */
	if (comm_rank == 0)
	{
		indices = new int [number_of_elements];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinCopyN(vecs_in[i]->getIndices(), numbers_of_elements[i], indices);
			indices += numbers_of_elements[i];
		}
		indices -= number_of_elements;
		MPI_Scatter(indices, number_of_elements, MPI_INT, NULL, 0, MPI_INT, 0, comm);
	}
	else
	{
		indices = new int [number_of_elements];
		MPI_Scatter(NULL, 0, MPI_INT, indices, number_of_elements, MPI_INT, 0, comm);
	}

	/** scatter values */
	if (comm_rank == 0)
	{
		values = new double [number_of_elements];
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinCopyN(vecs_in[i]->getElements(), vecs_in[i]->getNumElements(), values);
			values += vecs_in[i]->getNumElements();
		}
		values -= number_of_elements;
		MPI_Scatter(values, number_of_elements, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, comm);
	}
	else
	{
		values = new double [number_of_elements];
		MPI_Scatter(NULL, 0, MPI_DOUBLE, values, number_of_elements, MPI_DOUBLE, 0, comm);
	}

	/** re-construct output vectors */
	if (comm_rank > 0)
	{
		DSPdebugMessage("Rank %d: number_of_vectors %d\n", comm_rank, number_of_vectors);
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinPackedVector * vec = new CoinPackedVector(
					numbers_of_elements[i], indices, values);
			indices += numbers_of_elements[i];
			values  += numbers_of_elements[i];
			vecs_out.push_back(vec);
		}
		indices -= number_of_elements;
		values  -= number_of_elements;
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** MPI_Bcast function for vector<CoinPackedVectors*> */
DSP_RTN_CODE MPIbcastCoinPackedVectors(
		MPI_Comm comm,
		vector<CoinPackedVector*> & vecs)
{
#define FREE_MEMORY \
	FREE_ARRAY_PTR(numbers_of_vectors)  \
	FREE_ARRAY_PTR(numbers_of_elements) \
	FREE_ARRAY_PTR(indices)             \
	FREE_ARRAY_PTR(values)

	int   number_of_vectors;          /**< number of vectors in process */
	int * numbers_of_vectors = NULL;  /**< number of vectors per process */
	int   number_of_elements;         /**< number of elements in process */
	int * numbers_of_elements = NULL; /**< number of elements per vector in comm world */
	int * indices = NULL;
	double * values = NULL;

	BGN_TRY_CATCH

	int comm_rank;
	MPI_Comm_rank(comm, &comm_rank);

	/** broadcast the number of vectors */
	if (comm_rank == 0)
		number_of_vectors = vecs.size();
	MPI_Bcast(&number_of_vectors, 1, MPI_INT, 0, comm);

	/** broadcast the numbers of elements */
	numbers_of_elements = new int [number_of_vectors];
	if (comm_rank == 0)
	{
		for (int i = 0; i < number_of_vectors; ++i)
		{
			DSPdebugMessage("-> Process %d: vecs_in[%d]->getNumElements() %d\n",
					comm_rank, i, vecs[i]->getNumElements());
			numbers_of_elements[i] = vecs[i]->getNumElements();
		}
	}
	MPI_Bcast(numbers_of_elements, number_of_vectors, MPI_INT, 0, comm);

	/** get number of elements per process */
	number_of_elements = 0;
	for (int i = 0; i < number_of_vectors; ++i)
		number_of_elements += numbers_of_elements[i];

	/** broadcast indices */
	indices = new int [number_of_elements];
	if (comm_rank == 0)
	{
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinCopyN(vecs[i]->getIndices(), numbers_of_elements[i], indices);
			indices += numbers_of_elements[i];
		}
		indices -= number_of_elements;
	}
	MPI_Bcast(indices, number_of_elements, MPI_INT, 0, comm);

	/** broadcast values */
	values = new double [number_of_elements];
	if (comm_rank == 0)
	{
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinCopyN(vecs[i]->getElements(), vecs[i]->getNumElements(), values);
			values += vecs[i]->getNumElements();
		}
		values -= number_of_elements;
	}
	MPI_Bcast(values, number_of_elements, MPI_DOUBLE, 0, comm);

	/** re-construct output vectors */
	if (comm_rank > 0)
	{
		/** clear output vector */
		for (unsigned i = 0; i < vecs.size(); ++i)
			FREE_PTR(vecs[i]);
		vecs.clear();
		/** reserve memory */
		vecs.reserve(number_of_vectors);
		/** assign output vector */
		DSPdebugMessage("Rank %d: number_of_vectors %d\n", comm_rank, number_of_vectors);
		for (int i = 0; i < number_of_vectors; ++i)
		{
			CoinPackedVector * vec = new CoinPackedVector(
					numbers_of_elements[i], indices, values);
			indices += numbers_of_elements[i];
			values  += numbers_of_elements[i];
			vecs.push_back(vec);
		}
		indices -= number_of_elements;
		values  -= number_of_elements;
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** MPI_Scatter function for OsiCuts */
DSP_RTN_CODE MPIscatterOsiCuts(
		MPI_Comm comm,
		OsiCuts cuts_in,
		OsiCuts * cuts_out)
{
#define FREE_MEMORY     \
	FREE_ARRAY_PTR(lhs) \
	FREE_ARRAY_PTR(rhs)

	int comm_rank;
	int ncuts;
	vector<CoinPackedVector*> rows_in, rows_out;
	double * lhs = NULL;
	double * rhs = NULL;

	BGN_TRY_CATCH

	MPI_Comm_rank(comm, &comm_rank);

	if (comm_rank == 0)
	{
		/** allocate memory */
		ncuts = cuts_in.sizeCuts();
		lhs = new double [ncuts];
		rhs = new double [ncuts];

		/** parse cuts */
		for (int i = 0; i < ncuts; ++i)
		{
			OsiRowCut * rc = cuts_in.rowCutPtr(i);
			assert(rc);
			/** row vectors */
			rows_in.push_back(new CoinPackedVector(rc->row()));
			/** lhs */
			lhs[i] = rc->lb();
			/** rhs */
			rhs[i] = rc->ub();
		}
	}

	/** receive row vectors */
	MPIscatterCoinPackedVectors(comm, rows_in, rows_out);

	if (comm_rank > 0)
	{
		ncuts = rows_out.size();
		lhs = new double [ncuts];
		rhs = new double [ncuts];
	}


	/** scatter lhs */
	if (comm_rank == 0)
		MPI_Scatter(lhs, ncuts, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, comm);
	else
		MPI_Scatter(NULL, 0, MPI_DOUBLE, lhs, ncuts, MPI_DOUBLE, 0, comm);

	/** receive rhs */
	if (comm_rank == 0)
		MPI_Scatter(rhs, ncuts, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, comm);
	else
		MPI_Scatter(NULL, 0, MPI_DOUBLE, rhs, ncuts, MPI_DOUBLE, 0, comm);

	if (comm_rank > 0 && cuts_out != NULL)
	{
		/** construct cuts */
		for (int i = 0; i < ncuts; ++i)
		{
			OsiRowCut rc;
			rc.setRow(rows_out[i]);
			rc.setLb(lhs[i]);
			rc.setUb(rhs[i]);
			cuts_out->insert(rc);
		}
		/** free vector */
		for (unsigned i = 0; i < rows_out.size(); ++i)
			FREE_PTR(rows_out[i]);
	}

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}

/** MPI_Bcast function for OsiCuts */
DSP_RTN_CODE MPIbcastOsiCuts(
		MPI_Comm comm,
		OsiCuts * cuts)
{
#define FREE_MEMORY     \
	FREE_ARRAY_PTR(lhs) \
	FREE_ARRAY_PTR(rhs)

	int comm_rank;
	int ncuts;
	vector<CoinPackedVector*> rows;
	double * lhs = NULL;
	double * rhs = NULL;

	BGN_TRY_CATCH

	MPI_Comm_rank(comm, &comm_rank);

	if (comm_rank == 0)
	{
		/** allocate memory */
		ncuts = cuts->sizeCuts();
		lhs = new double [ncuts];
		rhs = new double [ncuts];

		/** parse cuts */
		for (int i = 0; i < ncuts; ++i)
		{
			OsiRowCut * rc = cuts->rowCutPtr(i);
			assert(rc);
			/** row vectors */
			rows.push_back(new CoinPackedVector(rc->row()));
			/** lhs */
			lhs[i] = rc->lb();
			/** rhs */
			rhs[i] = rc->ub();
		}
	}

	/** receive row vectors */
	MPIbcastCoinPackedVectors(comm, rows);
	if (comm_rank > 0)
	{
		ncuts = rows.size();
		lhs = new double [ncuts];
		rhs = new double [ncuts];
	}

	/** broadcast lhs */
	MPI_Bcast(lhs, ncuts, MPI_DOUBLE, 0, comm);
	/** broadcast rhs */
	MPI_Bcast(rhs, ncuts, MPI_DOUBLE, 0, comm);

	if (comm_rank > 0 && cuts != NULL)
	{
		/** clear output cuts */
		for (int i = 0; i < cuts->sizeCuts(); ++i)
		{
			OsiRowCut * rc = cuts->rowCutPtr(i);
			FREE_PTR(rc);
		}
		cuts->dumpCuts();
		/** construct cuts */
		for (int i = 0; i < ncuts; ++i)
		{
			OsiRowCut rc;
			rc.setRow(*rows[i]);
			rc.setLb(lhs[i]);
			rc.setUb(rhs[i]);
			cuts->insert(rc);
		}
	}
	/** free vector */
	for (unsigned i = 0; i < rows.size(); ++i)
		FREE_PTR(rows[i]);

	END_TRY_CATCH_RTN(FREE_MEMORY,DSP_RTN_ERR)

	FREE_MEMORY

	return DSP_RTN_OK;
#undef FREE_MEMORY
}
