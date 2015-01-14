include ../../make.inc

ZSRCS := $(wildcard ./*.f90)
ZOBJS := $(ZSRCS:.f90=.o)

all: srcs

srcs: $(ZOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I../../include

clean:
	@rm -f $(ZOBJS)	
