#!/bin/sh
#
# This script is used to run automated tests when "make check" is run.
#

PROG=$1

if [ -f $srcdir/$PROG.data ]
then
   ./$PROG < $srcdir/$PROG.data > $PROG.log
else
   ./$PROG > $PROG.log
fi
diff $PROG.log $srcdir/$PROG.output
