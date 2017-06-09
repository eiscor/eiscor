include ../../make.inc ../../.master.inc

DBMRKSRCS := $(wildcard *.f90)
DBMRKS := $(DBMRKSRCS:.f90=)

all: benchmarks

benchmarks:
	@$(MAKE) run

benchmarks_d: benchmarks

benchmarks_d_%: 
	$(eval DBMRKSRCS := $(patsubst benchmarks_%,benchmark_%,$@))
	$(eval DBMRKSRCS := $(wildcard $(DBMRKSRCS)*.f90))
	$(eval DBMRKS := $(DBMRKSRCS:.f90=))
	@$(MAKE) run DBMRKS="$(DBMRKS)"

run: $(DBMRKS)
	@$(foreach benchmark,$(DBMRKS),./$(benchmark) &&) echo 'End of double benchmarks!'

%:: %.f90
	$(FC) $(FFLAGS) -cpp $< -o $@ -I ../../include $(INSTALLDIR)/eiscor/lib/lib$(LIBNAME).$(SLIB).$(VERSION)  $(LIBS)

clean:
	@rm -f $(DBMRKS)	
