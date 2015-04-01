include ../../make.inc ../../.master.inc

LTESTSRCS := $(wildcard *.f90)
LTESTS := $(LTESTSRCS:.f90=)

all: tests

tests:
	@$(MAKE) run

tests_l: tests

tests_l_%: 
	$(eval LTESTSRCS := $(patsubst tests_%,test_%,$@))
	$(eval LTESTSRCS := $(wildcard $(LTESTSRCS)*.f90))
	$(eval LTESTS := $(LTESTSRCS:.f90=))
	@$(MAKE) run LTESTS="$(LTESTS)"

run: $(LTESTS)
	@$(foreach test,$(LTESTS),./$(test) &&) echo 'End of logical tests!'

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)  $(LIBS)

clean:
	@rm -f $(LTESTS)	
