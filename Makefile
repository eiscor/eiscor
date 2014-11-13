include make.inc

UTILDIR := ./src/utilities
ISRCDIR := ./src/integer
DSRCDIR := ./src/double
ZSRCDIR := ./src/complex_double

USRCS := $(wildcard $(UTILDIR)/*.f90)
UOBJS := $(patsubst $(UTILDIR)/%.f90,$(UTILDIR)/%.o,$(wildcard $(UTILDIR)/*.f90))

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

all: lib$(LIBNAME).so.$(VERSION)

examples: $(DEXS) $(ZEXS)
	$(DEXS) && $(ZEXS)

lib$(LIBNAME).so.$(VERSION): $(UOBJS) $(IOBJS) $(DOBJS) $(ZOBJS)
	$(FC) $(FFLAGS) -shared -o lib$(LIBNAME).so.$(VERSION) $(UOBJS) $(IOBJS) $(DOBJS) $(ZOBJS)
	
$(UOBJS): $(USRCS)
	make -C $(UTILDIR)
	
$(USRCS):

$(IOBJS): $(ISRCS)
	make -C $(ISRCDIR)
	
$(ISRCS):

$(DOBJS): $(DSRCS)
	make -C $(DSRCDIR)
	
$(DSRCS):

$(ZOBJS): $(ZSRCS)
	make -C $(ZSRCDIR)
	
$(ZSRCS):

$(DEXS): $(DEXSRCS) install
	make -C $(DEXDIR)
	
$(DEXSRCS):

$(ZEXS): $(ZEXSRCS) install
	make -C $(ZEXDIR)
	
$(ZEXSRCS):

install: lib$(LIBNAME).so.$(VERSION)
	mkdir -p $(INSTALLDIR)/$(LIBNAME)/lib &&\
	cp ./lib$(LIBNAME).so.$(VERSION) $(INSTALLDIR)/$(LIBNAME)/lib 

uninstall: clean
	rm -rf $(INSTALLDIR)/$(LIBNAME)

clean:
	make clean -C $(UTILDIR) &&\
	make clean -C $(ISRCDIR) &&\
	make clean -C $(DSRCDIR) &&\
	make clean -C $(ZSRCDIR) &&\
	make clean -C $(DEXDIR) &&\
	make clean -C $(ZEXDIR) &&\
	rm -f lib$(LIBNAME).so.$(VERSION)
	
	
