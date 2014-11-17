include ../../make.inc

DTESTSRCS := $(wildcard ./*.f90)
DTESTS := $(patsubst ./%.f90,./%,$(wildcard ./*.f90))

all: $(DTESTS)

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/$(LIBNAME)/lib/lib$(LIBNAME).so.$(VERSION) $(LIBS) && ./$@


clean:
	rm -f $(DTESTS)	
