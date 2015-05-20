include ../../make.inc

DSRCS := $(wildcard ./*.f90)
DOBJS := $(DSRCS:.f90=.o)

all: srcs

srcs: $(DOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I../../include

clean:
	@rm -f $(DOBJS)	
