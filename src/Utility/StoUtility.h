/*
 * StoUtility.h
 *
 *  Created on: Oct 7, 2014
 *      Author: kibaekkim
 */

#ifndef STOUTILITY_H_
#define STOUTILITY_H_

#include "mpi.h"

/** COIN */
#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"

using namespace std;

/** check whether solution is duplicate or not; return NULL if duplicate */
bool duplicateVector(
		CoinPackedVector * vec,
		vector<CoinPackedVector*> vecs);

/** MPI gather function for vector<CoinPackedVectors*> type. */
void MPIgatherCoinPackedVectors(
		MPI::Intracomm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out);

/** MPI scatter function for vector<CoinPackedVectors*> type. */
void MPIscatterCoinPackedVectors(
		MPI::Intracomm comm,
		vector<CoinPackedVector*> vecs_in,
		vector<CoinPackedVector*> & vecs_out);

/** MPI gather function for OsiCuts type. */
void MPIgatherOsiCuts(
		MPI::Intracomm comm,
		OsiCuts cuts_in,
		OsiCuts & cuts_out);

/** MPI scatter function for OsiCuts type. */
void MPIscatterOsiCuts(
		MPI::Intracomm comm,
		OsiCuts cuts_in,
		OsiCuts & cuts_out);

#endif /* STOUTILITY_H_ */
