include ../../make.inc

USRCS := $(wildcard ./*.f90)
UOBJS := $(patsubst ./%.f90,./%.o,$(wildcard ./*.f90))

all: $(UOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I../../include

clean:
	rm -f $(UOBJS)	
