include ../../make.inc

ZTESTSRCS := $(wildcard ./*.f90)
ZTESTS := $(patsubst ./%.f90,./%,$(wildcard ./*.f90))

all: $(ZTESTS)

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/$(LIBNAME)/lib/lib$(LIBNAME).so.$(VERSION) $(LIBS) && ./$@


clean:
	rm -f $(ZTESTS)	
