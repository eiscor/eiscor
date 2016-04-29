############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################

# directory where the library will be installed
INSTALLDIR := $(HOME)

# fortran compiler and flags
FC := gfortran
#FFLAGS := -O3 -std=f95
#FFLAGS := -pg
FFLAGS := -g -fbounds-check -Wall -fbacktrace -finit-real=nan

# libs required by compiler
LIBS := 




