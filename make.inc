############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################
LIBNAME := eiscor
MAJOR := 0
MINOR := 0
VERSION := $(MAJOR).$(MINOR)

# directory where the library will be installed
INSTALLDIR := $(HOME)

# fortran compiler and flags
FC := ifort 
FFLAGS := -std90 -heap-arrays

# libs required by compiler
LIBS := /opt/icc/composer_xe_2011_sp1.10.319/compiler/lib/intel64/libifport.so #-lifport






