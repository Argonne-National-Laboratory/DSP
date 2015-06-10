/*
 * StoMessage.h
 *
 *  Created on: Nov 26, 2014
 *      Author: kibaekkim
 */

#ifndef STOMESSAGE_H_
#define STOMESSAGE_H_

#include "CoinMessageHandler.hpp"

enum STO_Message
{
	STO_DD_END_GOOD,
	STO_DD_END_ITER,
	STO_DD_TIME,
	STO_DD_BEST_UB,
	STO_DD_BEST_LB,
	STO_DD_CURR_OBJVAL,
	STO_DD_TRY_UB,
	STO_DD_UPDATE_TR,
	STO_DD_ITER_INFO,
	STO_DD_SUBPROB_DETAIL,
	STO_DUMMY_END
};

class StoMessage: public CoinMessages
{
public:
	StoMessage(Language language = us_en);
};

#ifdef DSP_DEBUG

#define DSPdebug(x)        x
#define DSPdebugMessage    printf("[%s:%d] debug: ", __FILE__, __LINE__), printf

#else

#define DSPdebug(x)        while (false) x
#define DSPdebugMessage    while (false) printf

#endif

#endif /* STOMESSAGE_H_ */
