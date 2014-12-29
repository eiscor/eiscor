include ../make.inc

all: tests

tests:
	@$(MAKE) tests -C ./double
	@$(MAKE) tests -C ./complex_double
	
tests_d_%:
	@$(MAKE) $@ -C ./double
	
tests_z_%:
	@$(MAKE) $@ -C ./complex_double
	
clean:
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double
	
	
