include ../../make.inc ../../.master.inc

DEXSRCS := $(wildcard *.f90)
DEXS := $(DEXSRCS:.f90=)

all: examples

examples:
	@$(MAKE) run

examples_d: examples

examples_d_%: 
	$(eval DEXSRCS := $(patsubst examples_%,example_%,$@))
	$(eval DEXSRCS := $(wildcard $(DEXSRCS)*.f90))
	$(eval DEXS := $(DEXSRCS:.f90=))
	@$(MAKE) run DEXS="$(DEXS)"

run: $(DEXS)
	@$(foreach example,$(DEXS),./$(example) &&) echo 'End of double examples!'

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)  $(LIBS)

clean:
	@rm -f $(DEXS)	
