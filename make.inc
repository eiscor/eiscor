############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################

# directory where the library will be installed
INSTALLDIR := $(HOME)

# fortran compiler and flags
FC := gfortran
FFLAGS := -O3 -std=f95 -g -fbacktrace

# libs required by compiler
#LIBS := 
#LIBS := -llapack -lblas -lm 
#LIBS := /export/home1/NoCsBack/nalag/thomasm/lapack/lapack-3.5.0/lib/liblapack.a /export/home1/NoCsBack/nalag/thomasm/lapack/lapack-3.5.0/lib/libblas.a -lm
LIBS := ~/lapack/lapack-3.6.1/liblapack.a ~/lapack/BLAS-3.6.0/blas_LINUX.a -lm



