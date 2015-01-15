include ../../make.inc

LSRCS := $(wildcard ./*.f90)
LOBJS := $(LSRCS:.f90=.o)

all: srcs

srcs: $(LOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I../../include

clean:
	@rm -f $(LOBJS)	
