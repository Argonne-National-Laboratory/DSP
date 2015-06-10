/*
 * StoMessage.cpp
 *
 *  Created on: Nov 26, 2014
 *      Author: kibaekkim
 */

#include <cstring>
#include "Utility/StoMessage.h"

typedef struct
{
    STO_Message internalNumber;
    int externalNumber; // or continuation
    char detail;
    const char * message;
} Sto_message;

static Sto_message us_english[] =
{
		{STO_DD_END_GOOD,  0, 0, "Search completed - best upper bound %e, best lower bound %e, took %d iterations and solution time (cpu: %.2f wall: %.2f)"},
		{STO_DD_END_ITER,  1, 0, "Exiting on iteration limit - best upper bound %e, best lower bound %e, took %d iterations and solution time (cpu: %.2f wall: %.2f)"},
		{STO_DD_TIME,        10, 1, "Iteration %d, cpu (wall) time: master %.2f (%.2f)  subprob %.2f (%.2f)  other %.2f (%.2f)  total %.2f (%.2f)"},
		{STO_DD_BEST_UB,     11, 1, "Iteration %d, best upper bound %e optimality gap %.2f%%"},
		{STO_DD_BEST_LB,     12, 1, "Iteration %d, best lower bound %e approx lower bound %e cutoff gap %.2f%%"},
		{STO_DD_CURR_OBJVAL, 13, 1, "Iteration %d, Lagrangian dual value %e"},
		{STO_DD_TRY_UB,      14, 1, "Iteration %d, exploring upper bounds"},
		{STO_DD_UPDATE_TR,   15, 1, "Iteration %d, updating trust region with rho %e at cutoff %e"},
		{STO_DD_ITER_INFO,   16, 1, "Iteration %d, master %e dual bound %e gap %e"},
		{STO_DD_SUBPROB_DETAIL, 20, 2, "Scenario %3d, objval %e  current upper bound %e"},
		{STO_DUMMY_END, 999999, 0, ""}
};

StoMessage::StoMessage(Language language) :
	CoinMessages(sizeof(us_english) / sizeof(Sto_message))
{
	language_ = language;
	std::strcpy(source_, "Dsp");
	class_ = 1; /**< solver */
	Sto_message * message = us_english;

	while (message->internalNumber != STO_DUMMY_END)
	{
		CoinOneMessage oneMessage(message->externalNumber, message->detail,
				message->message);
		addMessage(message->internalNumber, oneMessage);
		message++;
	}
	// Put into compact form
	toCompact();

	// now override any language ones

	//switch (language) {

	//default:
	message = NULL;
	//  break;
	//}

	// replace if any found
	if (message)
	{
		while (message->internalNumber != STO_DUMMY_END)
		{
			replaceMessage(message->internalNumber, message->message);
			message++;
		}
	}
}

