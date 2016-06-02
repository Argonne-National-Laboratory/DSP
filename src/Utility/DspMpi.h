/*
 * DspUtility.h
 *
 *  Refactored on: April 8, 2016
 *     Created on: Oct 7, 2014
 *         Author: kibaekkim
 */

#ifndef DSPUTILITY_H_
#define DSPUTILITY_H_

#include "mpi.h"

#include "Utility/DspRtnCodes.h"

/** COIN */
#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"

using namespace std;

#define DSP_MPI_TAG_LB    0
#define DSP_MPI_TAG_CGBD  1
#define DSP_MPI_TAG_UB    2
#define DSP_MPI_TAG_SOLS  3

/** check whether solution is duplicate or not; return NULL if duplicate */
bool duplicateVector(
		CoinPackedVector * vec,
		vector<CoinPackedVector*> vecs);

/** MPI_Send for vector<CoinPackedVectors*> */
DSP_RTN_CODE MPIsendCoinPackedVectors(
		MPI::Intracomm comm,
		int to_rank,
		vector<CoinPackedVector*> vecs,
		int tag = MPI_ANY_TAG);

/** MPI_Recv for vector<CoinPackedVectors*> */
DSP_RTN_CODE MPIrecvCoinPackedVectors(
		MPI::Intracomm comm,
		int from_rank,
		vector<CoinPackedVector*> & vecs,
		int tag = MPI_ANY_TAG);

/** MPI gather function for vector<CoinPackedVectors*> type. */
DSP_RTN_CODE MPIgatherCoinPackedVectors(
		MPI::Intracomm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out);

/** MPI Allgather function for vector<CoinPackedVectors*> type. */
DSP_RTN_CODE MPIAllgatherCoinPackedVectors(
		MPI::Intracomm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out);

/** MPI_Scatter function for vector<CoinPackedVectors*> */
DSP_RTN_CODE MPIscatterCoinPackedVectors(
		MPI::Intracomm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out);

/** MPI_Send for OsiCuts */
DSP_RTN_CODE MPIsendOsiCuts(
		MPI::Intracomm comm,
		int to_rank,
		OsiCuts cuts,
		int tag = MPI_ANY_TAG);

/** MPI_Recv for OsiCuts */
DSP_RTN_CODE MPIrecvOsiCuts(
		MPI::Intracomm comm,
		int from_rank,
		OsiCuts & cuts,
		int tag = MPI_ANY_TAG);

/** MPI gather function for OsiCuts type. */
DSP_RTN_CODE MPIgatherOsiCuts(
		MPI::Intracomm comm,
		OsiCuts cuts_in,
		OsiCuts & cuts_out);

/** MPI Allgather function for OsiCuts type. */
DSP_RTN_CODE MPIAllgatherOsiCuts(
		MPI::Intracomm comm,
		OsiCuts cuts_in,
		OsiCuts & cuts_out);

/** MPI_Scatter function for OsiCuts */
DSP_RTN_CODE MPIscatterOsiCuts(
		MPI::Intracomm comm,
		OsiCuts cuts_in,
		OsiCuts * cuts_out);

#endif /* DSPUTILITY_H_ */
