# Build shared libraries for use with PARI/GP via sage
#
# Usage:
#   make              # build all .so files
#   make clean        # remove build artifacts
#
# Requires: sage (for PARI headers and library)

SAGE_LOCAL := $(shell sage -sh -c 'echo $$SAGE_LOCAL')

CC       = gcc
CFLAGS   = -O3 -march=native -Wall -fno-strict-aliasing -g -fPIC \
           -I"$(SAGE_LOCAL)/include"
LDFLAGS  = -shared -fPIC \
           -L$(SAGE_LOCAL)/lib -Wl,-rpath,$(SAGE_LOCAL)/lib
LIBS     = -lc -lm -lpari

TARGETS  = eqfminim.so orbmod2.so

all: $(TARGETS)

%.so: %.o
	$(CC) -o $@ $(LDFLAGS) $< $(LIBS)

%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $<

clean:
	rm -f *.o *.so

.PHONY: all clean
