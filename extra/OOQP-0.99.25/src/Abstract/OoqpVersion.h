/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef OOQPVERSIONH 
#define OOQPVERSIONH

#define OOQPVERSIONMAJOR 0
#define OOQPVERSIONMINOR 99
#define OOQPVERSIONPATCHLEVEL 25

#define OOQPVERSIONDATE "August 30, 2014"

#ifdef __cplusplus
extern "C"
{
#endif
  void printOoqpVersionString();
  void getOoqpVersionString( char buff[], int lbuff);
#ifdef __cplusplus
}
#endif

#endif
