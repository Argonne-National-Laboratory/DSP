/*
 * DspMacros.h
 *
 *  Refactored on: Apr 8, 2016
 *     Created on: Sep 23, 2014
 *         Author: kibaekkim
 */

#ifndef DSPMACROS_H_
#define DSPMACROS_H_

#include <cmath>

/*
 * Some print function
 */

#define PRINT_ARRAY(NUM,ARR)                                       \
	if (NUM > 0) {                                                 \
	std::cout << std::scientific;                                  \
	for (int __i = 0, __ii = 0; __i < NUM; ++__i) {                \
		if (fabs((ARR)[__i]) > 1.0e-10) {                          \
			std::cout << "  " << __i << " [" << (ARR)[__i] << "]"; \
			__ii++;                                                \
			if (__ii % 5 == 0) std::cout << std::endl;             \
		}                                                          \
	}                                                              \
	std::cout << std::fixed << std::endl;                          \
	}

#define PRINT_SPARSE_ARRAY(NUM,IND,VAL)                           \
	if (NUM > 0) {                                                \
	std::cout << std::scientific;                                 \
	for (int __i = 0; __i < NUM; ++__i) {                         \
		if (__i > 0 && __i % 5 == 0) std::cout << std::endl;      \
		std::cout << "  " << IND[__i] << " [" << VAL[__i] << "]"; \
	}                                                             \
	std::cout << std::fixed << std::endl;                         \
	}

#define PRINT_COIN_PACKED_VECTOR(VEC) \
	PRINT_SPARSE_ARRAY(VEC.getNumElements(), VEC.getIndices(), VEC.getElements())

#define PRINT_ARRAY_MSG(NUM,ARR,MSG)                    \
	std::cout << "Print <" << MSG << ">:" << std::endl; \
	PRINT_ARRAY(NUM,ARR)

#define PRINT_SPARSE_ARRAY_MSG(NUM,IND,VAL,MSG)         \
	std::cout << "Print <" << MSG << ">:" << std::endl; \
	PRINT_SPARSE_ARRAY(NUM,IND,VAL)

#define PRINT_COIN_PACKED_VECTOR_MSG(VEC,MSG)           \
	std::cout << "Print <" << MSG << ">:" << std::endl; \
	PRINT_COIN_PACKED_VECTOR(VEC)

#define PRINT_ARRAY_OBJECT(NUM,ARR,NAME)								\
	if (NUM > 0) {														\
	std::cout << std::scientific;                               		\
	for (int __i = 0, __ii = 0; __i < NUM; ++__i) {             		\
		if (fabs((ARR)[__i].NAME) > 1.0e-10) {                		\
		std::cout << "  " << __i << " [" << (ARR)[__i].NAME << "]";  	\
			__ii++;                                             		\
			if (__ii % 5 == 0) std::cout << std::endl;          		\
		}                                                       		\
	}                                                           		\
	std::cout << std::fixed << std::endl;                          		\
	}

#define PRINT_ARRAY_OBJECT_MSG(NUM,ARR,NAME,MSG)                    \
	std::cout << "Print <" << MSG << ">:" << std::endl; 			\
	PRINT_ARRAY_OBJECT(NUM,ARR,NAME)

/*
 * try-catch macros
 */

#define BGN_TRY_CATCH try { \
	CoinError::printErrors_ = true;

#define END_TRY_CATCH(STMT)                                                                \
	} catch (std::bad_alloc& ba) {                                                         \
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;                       \
		STMT                                                                               \
	} catch (const char * str) {                                                           \
		std::cerr << str << std::endl;                                                     \
		STMT                                                                               \
	} catch (std::exception & e) {                                                         \
		std::cerr << "Exception: " << e.what() << std::endl;                               \
		STMT                                                                               \
	} catch (...) {                                                                        \
		std::cerr << "Exception occurred at " << __FILE__ << ":" << __LINE__ << std::endl; \
		STMT                                                                               \
	}

#define END_TRY_CATCH_RTN(STMT,RTN_VAL)  \
	END_TRY_CATCH(STMT; return RTN_VAL;)


/*
 * Memory related macros
 */

#define FREE_PTR(PTR) \
	if (PTR) {        \
		delete PTR;   \
		PTR = NULL;   \
	}

#define FREE_2D_PTR(NUM,PTR)               \
	if (PTR) {                             \
		for (int _i = 0; _i < NUM; ++_i) { \
			FREE_PTR(PTR[_i])              \
		}                                  \
		delete [] PTR;                     \
		PTR = NULL;                        \
	}

#define FREE_ARRAY_PTR(PTR) \
	if (PTR) {              \
		delete [] PTR;      \
		PTR = NULL;         \
	}

#define FREE_2D_ARRAY_PTR(NUM,PTR)         \
	if (PTR) {                             \
		for (int _i = 0; _i < NUM; ++_i) { \
			FREE_ARRAY_PTR(PTR[_i])        \
		}                                  \
		delete [] PTR;                     \
		PTR = NULL;                        \
	}

/**
 * C Interface
 */

#define DSP_API_CHECK_ENV(RTN)                \
	if (env == NULL) {                        \
		printf("Error: Null API pointer.\n"); \
		return RTN;                           \
	}
#define DSP_API_CHECK_MODEL(RTN)                \
	DSP_API_CHECK_ENV(RTN)                      \
	if (env->model_ == NULL) {                    \
		printf("Error: Null model pointer.\n"); \
		return RTN;                             \
	}
#define DSP_API_CHECK_SOLVER(RTN)                \
	DSP_API_CHECK_ENV(RTN)                       \
	if (env->solver_ == NULL) {                  \
		printf("Error: Null solver pointer.\n"); \
		return RTN;                              \
	}

#endif /* DSPMACROS_H_ */
