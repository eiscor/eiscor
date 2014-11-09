############################################################
#
# ECORE (Eigensolvers based on unitary CORE transformations)
#
# Date:
#
#  October 2014
#
# Authors:
#
#  Jared L. Aurentz
#  Thomas Mach
#
# Related articles:
#
# 
#
############################################################
LIBNAME := ecore
MAJOR := 0
MINOR := 0
VERSION := $(MAJOR).$(MINOR)

INSTALLDIR := $(HOME)

FC := ifort 
FFLAGS := -std90 -heap-arrays

#LD_LIBRARY_PATH=/opt/icc/composer_xe_2011_sp1.10.319/compiler/lib/intel64:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH

LIBS := -lifport






