/*
 * DspSmps.c
 *
 *  Created on: Jun 23, 2016
 *      Author: kibaekkim
 */

#include <stdio.h>
#include <dlfcn.h>

#define ccall(lib, func_to_run, return_type, input_type, ...) \
	((return_type(*)input_type)dlsym(lib, func_to_run))(__VA_ARGS__)

int main(int argc, char * argv[])
{
	if (argc != 2)
	{
		printf("USAGE: DspSmps smps_file\n");
		goto error;
	}

	int rc = 0;
	const char * lib_file = "libDsp.so";
	const char * smps_file = argv[1];
	void * lib = NULL;
	void * dsp = NULL;

	lib = dlopen(lib_file, RTLD_LAZY);
	if (lib == NULL)
	{
		//printf("Error: Failed to open library <%s>\n", lib_file);
		fprintf(stderr, "dlopen failed: %s\n", dlerror());
		goto error;
	}

	/** create DSP environment */
	dsp = ccall(lib, "createEnv", void*, ());
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
	ccall(lib, "solveDd", void, (void*), dsp);

	/** finalize DSP */
	ccall(lib, "freeEnv", void, (void*), dsp);

	rc = dlclose(lib);
	if (rc != 0)
	{
		printf("Error: Failed to close library <%s>\n", lib_file);
		goto error;
	}

	return 0;

error:
	return 1;
}
