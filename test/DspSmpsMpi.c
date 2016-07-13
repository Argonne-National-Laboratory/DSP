/*
 * DspSmpsMpi.c
 *
 *  Created on: Jun 23, 2016
 *      Author: kibaekkim
 */

#include <stdio.h>
#include <dlfcn.h>
#include "mpi.h"

#define ccall(lib, func_to_run, return_type, input_type, ...) \
	((return_type(*)input_type)dlsym(lib, func_to_run))(__VA_ARGS__)

int main(int argc, char * argv[])
{
	if (argc != 2)
	{
		printf("USAGE: DspSmps smps_file\n");
		goto error;
	}

	/** initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm comm_ = MPI_COMM_WORLD;
	int comm_rank_;
	MPI_Comm_rank(comm_, &comm_rank_);

	int rc = 0;
	const char * lib_file = "libDsp.so";
	const char * smps_file = argv[1];

	void * lib = dlopen(lib_file, RTLD_LAZY);
	if (lib == NULL)
	{
		printf("Error: Failed to open library <%s>\n", lib_file);
		goto error;
	}

	/** create DSP environment */
	void * dsp = ccall(lib, "createEnv", void*, ());
	if (dsp == NULL)
	{
		printf("Error: Failed to create DSP environment\n");
		goto error;
	}

	/** read SMPS file */
	ccall(lib, "readSmps", void, (void*,const char*), dsp, smps_file);

	/** set parameters */
	ccall(lib, "readParamFile", void, (void*,const char*), dsp, "params.txt");

	/** solve */
	ccall(lib, "solveDdMpi", void, (void*,MPI_Comm), dsp, comm);

	/** finalize DSP */
	ccall(lib, "freeEnv", void, (void*), dsp);

	rc = dlclose(lib);
	if (rc != 0)
	{
		printf("Error: Failed to close library <%s>\n", lib_file);
		goto error;
	}

	/** finalize MPI */
	MPI_Finalize();

	return 0;

error:
	return 1;
}
