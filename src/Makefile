include ../make.inc

all: srcs

srcs:
	@$(MAKE) -C ./utilities
	@$(MAKE) -C ./integer
	@$(MAKE) -C ./double
	@$(MAKE) -C ./complex_double

clean:
	@$(MAKE) clean -C ./utilities
	@$(MAKE) clean -C ./integer
	@$(MAKE) clean -C ./double
	@$(MAKE) clean -C ./complex_double

