############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################

# directory where the library will be installed
PREFIX := $(HOME)

# fortran compiler and flags
FC := gfortran
FFLAGS := -O3 -std=f95

# libs required by compiler
LIBS := 

# Installation prefix and library directory
PREFIX ?= /usr/local
LIBDIR ?= $(PREFIX)/lib

# Detect OS to name shared library accordingly
ifeq ($(OS),Windows_NT)
    IS_WINDOWS := 1
else
    UNAME_S := $(shell uname -s)
    ifeq ($(filter MSYS_NT% MINGW% CYGWIN_NT%,$(UNAME_S)),)
        IS_WINDOWS := 0
    else
        IS_WINDOWS := 1
    endif
endif

ifeq ($(IS_WINDOWS),1)
    SHARED_LIB = $(LIBNAME)-$(VERSION).dll
    INSTALL_LIB = cp $(SHARED_LIB) $(LIBDIR)/$(SHARED_LIB)
else
    ifeq ($(UNAME_S),Linux)
        SHARED_LIB_REAL = lib$(LIBNAME).so.$(VERSION)
        SHARED_LIB = lib$(LIBNAME).so
        INSTALL_LIB = install -m 755 $(SHARED_LIB) $(LIBDIR)/$(SHARED_LIB_REAL); \
            ln -sf $(LIBDIR)/$(SHARED_LIB_REAL) $(LIBDIR)/$(SHARED_LIB); \
            ln -sf $(LIBDIR)/$(SHARED_LIB_REAL) $(LIBDIR)/lib$(LIB_NAME).so.$(MAJOR); \
        patchelf --set-rpath $(LIBDIR) $(LIBDIR)/$(SHARED_LIB_REAL)
    else ifeq ($(UNAME_S),Darwin)
        SHARED_LIB_REAL = lib$(LIBNAME).$(VERSION).dylib
        SHARED_LIB = lib$(LIBNAME).dylib
        INSTALL_LIB = install -m 755 $(SHARED_LIB) $(LIBDIR)/$(SHARED_LIB_REAL); \
            ln -sf $(LIBDIR)/$(SHARED_LIB_REAL) $(LIBDIR)/$(SHARED_LIB); \
            install_name_tool -id $(LIBDIR)/$(SHARED_LIB_REAL) $(LIBDIR)/$(SHARED_LIB_REAL)
    endif
endif


