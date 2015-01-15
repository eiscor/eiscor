include ../make.inc

all: srcs

srcs:
	@$(MAKE) $@ -C ./utilities
	@$(MAKE) $@ -C ./logical
	@$(MAKE) $@ -C ./double
	@$(MAKE) $@ -C ./complex_double

clean:
	@$(MAKE) clean -C ./utilities
	@$(MAKE) clean -C ./logical
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double

