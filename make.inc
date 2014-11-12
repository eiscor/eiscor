############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################
LIBNAME := eiscor
MAJOR := 0
MINOR := 0
VERSION := $(MAJOR).$(MINOR)

INSTALLDIR := $(HOME)

FC := ifort 
FFLAGS := -std90 -heap-arrays

LIBS := -lifport






