include ../make.inc ../.master.inc

all: tests

tests:
	@$(MAKE) $@ -C ./double
	@$(MAKE) $@ -C ./complex_double
	
tests_d:
	@$(MAKE) $@ -C ./double
	
tests_z:
	@$(MAKE) $@ -C ./complex_double
	
tests_d_%:
	@$(MAKE) $@ -C ./double
	
tests_z_%:
	@$(MAKE) $@ -C ./complex_double
	
clean:
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double
	
	
