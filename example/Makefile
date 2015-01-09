include ../make.inc ../.master.inc

all: examples

examples:
	@$(MAKE) $@ -C ./double
	@$(MAKE) $@ -C ./complex_double

examples_d:
	@$(MAKE) $@ -C ./double

examples_z:
	@$(MAKE) $@ -C ./complex_double

examples_d_%:
	@$(MAKE) $@ -C ./double

examples_z_%:
	@$(MAKE) $@ -C ./complex_double

clean:
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double

