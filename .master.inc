# current version of eiscor
LIBNAME := eiscor
MAJOR := 0
MINOR := 1
PATCH := 0
VERSION := $(MAJOR).$(MINOR).$(PATCH)

# change library extension based on OS
ifeq ($(OS), Windows_NT)
	SLIB = dll
else
	UNAME := $(shell uname)
	ifeq ($(UNAME), Darwin)
		SLIB = dylib
	else
		SLIB = so
	endif
endif
