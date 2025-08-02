include make.inc .master.inc

SRCS := $(wildcard ./src/*/*.f90)
OBJS := $(SRCS:.f90=.o)

all: $(SHARED_LIB)

install: $(SHARED_LIB)
	@mkdir -p $(LIBDIR)
	$(INSTALL_LIB)

example%:
	@$(MAKE) $@ -C ./example

test%:
	@$(MAKE) $@ -C ./test

benchmark%:
	@$(MAKE) $@ -C ./benchmark

$(SHARED_LIB): srcs
	$(FC) $(FFLAGS) -shared -o $(SHARED_LIB) $(OBJS)

srcs:
	@$(MAKE) $@ -C ./src

uninstall: clean
	@rm -rf $(INSTALL_LIB)

clean:
	@$(MAKE) clean -C ./src
	@$(MAKE) clean -C ./example
	@$(MAKE) clean -C ./test
	@$(MAKE) clean -C ./benchmark
	@rm -f $(SHARED_LIB)


