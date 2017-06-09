include ../make.inc ../.master.inc

all: tests

tests:
	@$(MAKE) $@ -C ./logical
	@$(MAKE) $@ -C ./double
	@$(MAKE) $@ -C ./complex_double

tests_l:
	@$(MAKE) $@ -C ./logical

tests_d:
	@$(MAKE) $@ -C ./double

tests_z:
	@$(MAKE) $@ -C ./complex_double

tests_l_%:
	@$(MAKE) $@ -C ./logical

tests_d_%:
	@$(MAKE) $@ -C ./double

tests_z_%:
	@$(MAKE) $@ -C ./complex_double

clean:
	@$(MAKE) clean -C ./logical
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double
