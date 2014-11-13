include ../../make.inc

ISRCS := $(wildcard ./*.f90)
IOBJS := $(patsubst ./%.f90,./%.o,$(wildcard ./*.f90))

all: $(IOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I ../../include

clean:
	rm -f $(IOBJS)	
