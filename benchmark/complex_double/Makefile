include ../../make.inc ../../.master.inc

ZBMRKSRCS := $(wildcard *.f90)
ZBMRKS := $(ZBMRKSRCS:.f90=)

all: benchmarks

benchmarks:
	@$(MAKE) run

benchmarks_z: benchmarks

benchmarks_z_%: 
	$(eval ZBMRKSRCS := $(patsubst benchmarks_%,benchmark_%,$@))
	$(eval ZBMRKSRCS := $(wildcard $(ZBMRKSRCS)*.f90))
	$(eval ZBMRKS := $(ZBMRKSRCS:.f90=))
	@$(MAKE) run ZBMRKS="$(ZBMRKS)"

run: $(ZBMRKS)
	@$(foreach benchmark,$(ZBMRKS),./$(benchmark) &&) echo 'End of complex_double benchmarks!'

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION) $(LIBS)

clean:
	@rm -f $(ZBMRKS)	
