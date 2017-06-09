include ../../make.inc ../../.master.inc

ZEXSRCS := $(wildcard *.f90)
ZEXS := $(ZEXSRCS:.f90=)

all: examples

examples:
	@$(MAKE) run

examples_z: examples

examples_z_%: 
	$(eval ZEXSRCS := $(patsubst examples_%,example_%,$@))
	$(eval ZEXSRCS := $(wildcard $(ZEXSRCS)*.f90))
	$(eval ZEXS := $(ZEXSRCS:.f90=))
	@$(MAKE) run ZEXS="$(ZEXS)"

run: $(ZEXS)
	@$(foreach example,$(ZEXS),./$(example) &&) echo 'End of complex_double examples!'

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)  $(LIBS)

clean:
	@rm -f $(ZEXS)	
