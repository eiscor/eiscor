include ../../make.inc ../../.master.inc

DTESTSRCS := $(wildcard *.f90)
DTESTS := $(DTESTSRCS:.f90=)

all: tests

tests:
	@$(MAKE) run

tests_d: tests

tests_d_%: 
	$(eval DTESTSRCS := $(patsubst tests_%,test_%,$@))
	$(eval DTESTSRCS := $(wildcard $(DTESTSRCS)*.f90))
	$(eval DTESTS := $(DTESTSRCS:.f90=))
	@$(MAKE) run DTESTS="$(DTESTS)"

run: $(DTESTS)
	@$(foreach test,$(DTESTS),./$(test) &&) echo 'End of double tests!'

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)  $(LIBS)

clean:
	@rm -f $(DTESTS)	
