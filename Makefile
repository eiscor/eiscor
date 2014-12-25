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

USRCDIR := ./src/utilities
ISRCDIR := ./src/integer
DSRCDIR := ./src/double
ZSRCDIR := ./src/complex_double

USRCS := $(wildcard $(USRCDIR)/*.f90)
UOBJS := $(patsubst $(USRCDIR)/%.f90,$(USRCDIR)/%.o,$(wildcard $(USRCDIR)/*.f90))

ISRCS := $(wildcard $(ISRCDIR)/*.f90)
IOBJS := $(patsubst $(ISRCDIR)/%.f90,$(ISRCDIR)/%.o,$(wildcard $(ISRCDIR)/*.f90))

DSRCS := $(wildcard $(DSRCDIR)/*.f90)
DOBJS := $(patsubst $(DSRCDIR)/%.f90,$(DSRCDIR)/%.o,$(wildcard $(DSRCDIR)/*.f90))

ZSRCS := $(wildcard $(ZSRCDIR)/*.f90)
ZOBJS := $(patsubst $(ZSRCDIR)/%.f90,$(ZSRCDIR)/%.o,$(wildcard $(ZSRCDIR)/*.f90))

DEXDIR := ./examples/double
ZEXDIR := ./examples/complex_double

DEXSRCS := $(wildcard $(DEXDIR)/*.f90)
DEXS := $(patsubst $(DEXDIR)/%.f90,$(DEXDIR)/%,$(wildcard $(DEXDIR)/*.f90))

ZEXSRCS := $(wildcard $(ZEXDIR)/*.f90)
ZEXS := $(patsubst $(ZEXDIR)/%.f90,$(ZEXDIR)/%,$(wildcard $(ZEXDIR)/*.f90))

DTESTDIR := ./tests/double
ZTESTDIR := ./tests/complex_double

DTESTSRCS := $(wildcard $(DTESTDIR)/*.f90)
DTESTS := $(patsubst $(DTESTDIR)/%.f90,$(DTESTDIR)/%,$(wildcard $(DTESTDIR)/*.f90))

ZTESTSRCS := $(wildcard $(ZTESTDIR)/*.f90)
ZTESTS := $(patsubst $(ZTESTDIR)/%.f90,$(ZTESTDIR)/%,$(wildcard $(ZTESTDIR)/*.f90))

all: lib$(LIBNAME).$(SLIB).$(VERSION)

install: lib$(LIBNAME).$(SLIB).$(VERSION)
	mkdir -p $(INSTALLDIR)/$(LIBNAME)/lib &&\
	cp ./lib$(LIBNAME).$(SLIB).$(VERSION) $(INSTALLDIR)/$(LIBNAME)/lib

examples: $(DEXS) $(ZEXS)
	$(foreach ex,$(DEXS),$(ex) &&) \
	$(foreach ex,$(ZEXS),$(ex) &&) \
	echo 'End of examples!'

tests: $(DTESTS) $(ZTESTS)
	$(foreach test,$(DTESTS),./tests/double/$(test) &&) \
	$(foreach test,$(ZTESTS),./tests/complex_double/$(test) &&) \
	echo 'End of tests!'

lib$(LIBNAME).$(SLIB).$(VERSION): $(UOBJS) $(IOBJS) $(DOBJS) $(ZOBJS)
	$(FC) $(FFLAGS) -shared -o lib$(LIBNAME).$(SLIB).$(VERSION) $(UOBJS) $(IOBJS) $(DOBJS) $(ZOBJS) 

$(UOBJS) $(IOBJS) $(DOBJS) $(ZOBJS): $(USRCS) $(ISRCS) $(DSRCS) $(ZSRCS)
	make -C $(USRCDIR) &&\
	make -C $(ISRCDIR) &&\
	make -C $(DSRCDIR) &&\
	make -C $(ZSRCDIR)

$(USRCS) $(ISRCS) $(DSRCS) $(ZSRCS):

$(DEXS) $(ZEXS): $(DEXSRCS) $(ZEXSRCS) lib$(LIBNAME).$(SLIB).$(VERSION)
	make -C $(DEXDIR) &&\
	make -C $(ZEXDIR)

$(DEXSRCS) $(ZEXSRCS):

$(DTESTS) $(ZTESTS): $(DTESTSRCS) $(ZTESTSRCS) lib$(LIBNAME).$(SLIB).$(VERSION)
	make -C $(DTESTDIR) &&\
	make -C $(ZTESTDIR)

$(DTESTSRCS) $(ZTESTSRCS):

uninstall: clean
	rm -rf $(INSTALLDIR)/$(LIBNAME)

clean:
	make clean -C $(USRCDIR) &&\
	make clean -C $(ISRCDIR) &&\
	make clean -C $(DSRCDIR) &&\
	make clean -C $(ZSRCDIR) &&\
	make clean -C $(DEXDIR) &&\
	make clean -C $(ZEXDIR) &&\
	make clean -C $(DTESTDIR) &&\
	make clean -C $(ZTESTDIR) &&\
	rm -f lib$(LIBNAME).so.$(VERSION)


