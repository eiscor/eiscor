include ../../make.inc

ZSRCS := $(wildcard ./*.f90)
ZOBJS := $(patsubst ./%.f90,./%.o,$(wildcard ./*.f90))

all: $(ZOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I../../include

clean:
	rm -f $(ZOBJS)	
