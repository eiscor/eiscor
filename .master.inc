# current version of eiscor
LIBNAME := eiscor
MAJOR := 0
MINOR := 2
PATCH := 0
VERSION := $(MAJOR).$(MINOR).$(PATCH)

# change library extension based on OS
UNAME := $(shell uname)
ifeq ($(findstring Windows_NT,$(UNAME)),Windows_NT)
	SLIB := dll
else ifeq ($(findstring Darwin,$(UNAME)),Darwin)
	SLIB := dylib
else ifeq ($(findstring darwin,$(UNAME)),darwin)
	SLIB := dylib
else
	SLIB := so
endif
