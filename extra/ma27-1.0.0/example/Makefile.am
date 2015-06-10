AM_FCFLAGS = -I $(top_builddir)/src

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	ma27ss \
	ma27ds

check_PROGRAMS = \
	ma27ss \
	ma27ds

ma27ss_SOURCES = ma27ss.f
ma27ss_LDADD = $(top_builddir)/src/libma27.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma27ss.output ma27ss.data
ma27ds_SOURCES = ma27ds.f
ma27ds_LDADD = $(top_builddir)/src/libma27.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += ma27ds.output ma27ds.data
