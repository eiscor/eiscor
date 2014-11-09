include ../../make.inc

DSRCS := $(wildcard ./*.f90)
DOBJS := $(patsubst ./%.f90,./%.o,$(wildcard ./*.f90))

all: $(DOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I../../include

clean:
	rm -f $(DOBJS)	
