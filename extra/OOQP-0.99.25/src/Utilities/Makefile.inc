UTILDIR     = $(srcdir)/Utilities

UTILSOURCES = $(UTILDIR)/storage.C $(UTILDIR)/dsytrs.f $(UTILDIR)/dsytrf.f \
    $(UTILDIR)/dgetrs.f $(UTILDIR)/dgetrf.f \
    $(UTILDIR)/dpotrs.f	$(UTILDIR)/dpotrf.f \
    $(UTILDIR)/drand.C $(UTILDIR)/sort.C

UTILOBJECTS = $(UTILDIR)/dsytrs.o $(UTILDIR)/dsytrf.o \
	$(UTILDIR)/dgetrs.o $(UTILDIR)/dgetrf.o \
	$(UTILDIR)/dpotrs.o $(UTILDIR)/dpotrf.o	\
	$(UTILDIR)/drand.o $(UTILDIR)/sort.o

UTILHEADERS = $(patsubst $(UTILDIR)/%,$(localincludedir)/%,\
       $(wildcard $(UTILDIR)/*.h))