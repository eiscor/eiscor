include ../../make.inc ../../.master.inc

ZTESTSRCS := $(wildcard *.f90)
ZTESTS := $(ZTESTSRCS:.f90=)

all: tests

tests:
	@$(MAKE) run

tests_z: tests

tests_z_%: 
	$(eval ZTESTSRCS := $(patsubst tests_%,test_%,$@))
	$(eval ZTESTSRCS := $(wildcard $(ZTESTSRCS)*.f90))
	$(eval ZTESTS := $(ZTESTSRCS:.f90=))
	@$(MAKE) run ZTESTS="$(ZTESTS)"

run: $(ZTESTS)
	@$(foreach test,$(ZTESTS),./$(test) &&) echo 'End of complex_double tests!'

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION) $(LIBS)

clean:
	@rm -f $(ZTESTS)	
