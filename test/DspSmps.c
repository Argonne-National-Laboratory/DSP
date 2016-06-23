/*
 * DspSmps.c
 *
 *  Created on: Jun 23, 2016
 *      Author: kibaekkim
 */

#include <stdio.h>
#include "dbg.h"
#include <dlfcn.h>
#include "mpi.h"

#define ccall(lib, func_to_run, return_type, input_type, ...) \
	((return_type(*)input_type)dlsym(lib, func_to_run))(__VA_ARGS__)

int main(int argc, char * argv[])
{
	if (argc != 3)
	{
		printf("USAGE: DspSmps smps_file output_file\n");
		goto error;
	}

	/** initialize MPI */
	MPI_Init(&argc, &argv);

	int rc = 0;
	const char * lib_file = "libDsp.so";
	const char * smps_file = argv[1];
	const char * param_file = argv[2];

	void * lib = dlopen(lib_file, RTLD_LAZY);
	check(lib != NULL, "Failed to open the library %s: %s", lib_file, dlerror());

	/** create DSP environment */
	void * dsp = ccall(lib, "createEnv", void*, ());

	/** read SMPS file */
	ccall(lib, "readSmps", void, (void*,const char*), dsp, smps_file);

	/** set parameters */
	ccall(lib, "readParamFile", void, (void*,const char*), dsp, param_file);

	/** solve */
	ccall(lib, "solveDd", void, (void*,MPI_Comm), dsp, MPI_COMM_WORLD);

	/** finalize DSP */
	ccall(lib, "freeEnv", void, (void*), dsp);

	rc = dlclose(lib);
	check(rc == 0, "Failed to close %s", lib_file);

	/** finalize MPI */
	MPI_Finalize();

	return 0;

error:
	return 1;
}
