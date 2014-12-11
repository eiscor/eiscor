include ../../make.inc

DEXSRCS := $(wildcard ./*.f90)
DEXS := $(patsubst ./%.f90,./%,$(wildcard ./*.f90))

all: $(DEXS)

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(EISCORDIR)/lib$(LIBNAME).so.$(VERSION) $(LIBS)


clean:
	rm -f $(DEXS)	
