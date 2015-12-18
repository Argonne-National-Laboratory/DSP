/*
 * StoMessage.h
 *
 *  Created on: Nov 26, 2014
 *      Author: kibaekkim
 */

#ifndef STOMESSAGE_H_
#define STOMESSAGE_H_

#include <stdarg.h>
#include <stdio.h>

class StoMessage
{
public:
	StoMessage(int logLevel): logLevel_(logLevel) {}

	void print(int level, const char *fmt, ...)
	{
		if (level <= logLevel_)
		{
			va_list args;
			va_start(args, fmt);
			vfprintf(stderr, fmt, args);
			va_end(args);
		}
	}

	int logLevel_;
};

#ifdef DSP_DEBUG

#define DSPdebug(x)        x
#define DSPdebugMessage    printf("[%s:%d] debug: ", __FILE__, __LINE__), printf

#else

#define DSPdebug(x)        while (false) x
#define DSPdebugMessage    while (false) printf

#endif

#endif /* STOMESSAGE_H_ */
