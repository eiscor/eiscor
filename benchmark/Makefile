include ../make.inc ../.master.inc

all: benchmarks

benchmarks:
	@$(MAKE) $@ -C ./double
	@$(MAKE) $@ -C ./complex_double

benchmarks_d:
	@$(MAKE) $@ -C ./double

benchmarks_z:
	@$(MAKE) $@ -C ./complex_double

benchmarks_d_%:
	@$(MAKE) $@ -C ./double

benchmarks_z_%:
	@$(MAKE) $@ -C ./complex_double

clean:
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double
