dnl Process this file with autoconf to produce a configure script.
AC_INIT
AC_CONFIG_SRCDIR([daxpy.f])

dnl Checks for programs.

dnl If FC is set, but F77 is not, assume that FC points to the F77
dnl compiler.
if test -z "$F77"; then
   if test -n "$FC"; then
      F77="$FC"
      export F77
   fi
fi

dnl Use native F77 first
if test -z "$F77"; then
   AC_CHECK_PROG(F77,f77,f77)
fi
dnl And optimize, unless told otherwise
if test -z "$FFLAGS"; then FFLAGS=-O; fi
AC_PROG_F77

dnl Checks for libraries.

dnl Checks for header files.

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
