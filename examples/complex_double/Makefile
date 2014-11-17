include ../../make.inc

ZEXSRCS := $(wildcard ./*.f90)
ZEXS := $(patsubst ./%.f90,./%,$(wildcard ./*.f90))

all: $(ZEXS)

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/$(LIBNAME)/lib/lib$(LIBNAME).so.$(VERSION) $(LIBS) && ./$@

clean:
	rm -f $(ZEXS)	
