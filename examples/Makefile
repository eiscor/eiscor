include ../make.inc

all: examples

examples:
	@$(MAKE) examples -C ./double
	@$(MAKE) examples -C ./complex_double
	
examples_d_%:
	@$(MAKE) $@ -C ./double
	
examples_z_%:
	@$(MAKE) $@ -C ./complex_double
	
clean:
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double
	
	
