include ../../make.inc

ZEXSRCS := $(wildcard ./*.f90)
ZEXS := $(patsubst ./%.f90,./%,$(wildcard ./*.f90))

all: $(ZEXS)

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ $(INSTALLDIR)/$(LIBNAME)/lib/lib$(LIBNAME).so $(LIBS)

clean:
	rm -f $(ZEXS)	
