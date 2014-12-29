include make.inc

ifeq ($(OS), Windows_NT)
	SLIB = dll
else
	UNAME := $(shell uname)
	ifeq ($(UNAME), Darwin)
		SLIB = dylib
	else
		SLIB = so
	endif
endif

SRCS := $(wildcard ./src/*/*.f90)
OBJS := $(SRCS:.f90=.o)

all: lib$(LIBNAME).$(SLIB).$(VERSION)

install: lib$(LIBNAME).$(SLIB).$(VERSION)
	@mkdir -p $(INSTALLDIR)/$(LIBNAME)/lib
	@cp ./lib$(LIBNAME).$(SLIB).$(VERSION) $(INSTALLDIR)/$(LIBNAME)/lib

example%: install
	@$(MAKE) $@ -C ./examples

test%: install
	@$(MAKE) $@ -C ./tests

lib$(LIBNAME).$(SLIB).$(VERSION): $(OBJS)
	$(FC) $(FFLAGS) -shared -o lib$(LIBNAME).$(SLIB).$(VERSION) $(OBJS) 

$(OBJS): $(SRCS)
	@$(MAKE) -C ./src

$(SRCS):

uninstall: clean
	@rm -rf $(INSTALLDIR)/$(LIBNAME)

clean:
	@$(MAKE) clean -C ./src
	@$(MAKE) clean -C ./examples
	@$(MAKE) clean -C ./tests
	@rm -f lib$(LIBNAME).so.$(VERSION)


