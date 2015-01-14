include ../../make.inc

USRCS := $(wildcard ./*.f90)
UOBJS := $(USRCS:.f90=.o)

all: srcs

srcs: $(UOBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -cpp -fPIC -c $< -I../../include

clean:
	@rm -f $(UOBJS)	
