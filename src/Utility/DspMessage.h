/*
 * DspMessage.h
 *
 *  Refactored on: Apr 8, 2016
 *     Created on: Nov 26, 2014
 *         Author: kibaekkim
 */

#ifndef DSPMESSAGE_H_
#define DSPMESSAGE_H_

#include <stdarg.h>
#include <stdio.h>
#include <cmath>
#include "CoinPackedVector.hpp"

class DspMessage
{
public:
	DspMessage(int logLevel): logLevel_(logLevel)
	{
		setbuf(stdout, NULL);
	}

	void print(int level, const char *fmt, ...)
	{
		if (level <= logLevel_)
		{
			va_list args;
			va_start(args, fmt);
			vfprintf(stdout, fmt, args);
			va_end(args);
		}
	}

	int logLevel_;

	static void printArray(const CoinPackedVector * values)
	{
		for (int i = 0, j = 0; i < values->getNumElements(); ++i)
		{
			if (j > 0 && j % 5 == 0) printf("\n");
			printf("  [%6d] %+e", values->getIndices()[i], values->getElements()[i]);
			j++;
		}
		printf("\n");
	}

	static void printArray(int n, const int * indices, const double * values)
	{
		for (int i = 0, j = 0; i < n; ++i)
		{
			if (j > 0 && j % 5 == 0) printf("\n");
			printf("  [%6d] %+e", indices[i], values[i]);
			j++;
		}
		printf("\n");
	}

	static void printArray(int n, const double * values)
	{
		for (int i = 0, j = 0; i < n; ++i)
		{
			if (fabs(values[i]) < 1.0e-10) continue;
			if (j > 0 && j % 5 == 0) printf("\n");
			printf("  [%6d] %+e", i, values[i]);
			j++;
		}
		printf("\n");
	}

	static void printArray(int n, const int * values)
	{
		for (int i = 0, j = 0; i < n; ++i)
		{
			if (fabs(values[i]) < 1.0e-10) continue;
			if (j > 0 && j % 5 == 0) printf("\n");
			printf("  [%6d] %+10d", i, values[i]);
			j++;
		}
		printf("\n");
	}
};

#ifdef DSP_DEBUG

#define DSPdebug(x)        x
#define DSPdebugMessage    printf("[%s:%d] debug: ", __FILE__, __LINE__), printf

#else/* DSP_DEBUG */

#define DSPdebug(x)        while (false) x
#define DSPdebugMessage    while (false) printf

#endif/* DSP_DEBUG */

#ifdef DSP_DEBUG2

#define DSPdebug2(x)       x
#define DSPdebugMessage2   printf("[%s:%d] debug: ", __FILE__, __LINE__), printf

#else/* DSP_DEBUG2 */

#define DSPdebug2(x)       while (false) x
#define DSPdebugMessage2   while (false) printf

#endif/* DSP_DEBUG2 */

#endif /* DSPMESSAGE_H_ */
