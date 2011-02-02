DEFORMATION_MPIR_LIB_DIR=/home/suser/MPIR/mpir-2.1.1//.libs
DEFORMATION_MPIR_INCLUDE_DIR=/home/suser/MPIR/mpir-2.1.1/
DEFORMATION_MPFR_LIB_DIR=/home/suser/MPFR/mpfr-3.0.0//.libs
DEFORMATION_MPFR_INCLUDE_DIR=/home/suser/MPFR/mpfr-3.0.0/
DEFORMATION_FLINT_LIB_DIR=/home/suser/FLINT/flint2-seb
DEFORMATION_FLINT_INCLUDE_DIR=/home/suser/FLINT/flint2-seb

DEFORMATION_LIB=libdeformation.so
CC=gcc
CFLAGS=-O2 -g -ansi -pedantic -Wall -funroll-loops -Wno-unused
PREFIX=~/DEFORMATION/

LIBS=-L$(CURDIR) -L$(DEFORMATION_MPIR_LIB_DIR) -L$(DEFORMATION_MPFR_LIB_DIR) -L$(DEFORMATION_FLINT_LIB_DIR) -ldeformation -lflint -lmpir -lmpfr -lm
LIBS2=-L$(DEFORMATION_MPIR_LIB_DIR) -L$(DEFORMATION_MPFR_LIB_DIR) -L$(DEFORMATION_FLINT_LIB_DIR) -lflint -lmpir -lmpfr -lm
INCS=-I$(CURDIR) -I$(DEFORMATION_MPIR_INCLUDE_DIR) -I$(DEFORMATION_MPFR_INCLUDE_DIR) -I$(DEFORMATION_FLINT_INCLUDE_DIR)

LD_LIBRARY_PATH:=${CURDIR}:${DEFORMATION_FLINT_LIB_DIR}:${DEFORMATION_MPFR_LIB_DIR}:${DEFORMATION_MPIR_LIB_DIR}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
export

SOURCES = $(wildcard *.c)

HEADERS = $(wildcard *.h)

OBJS = $(patsubst %.c, %.o, $(SOURCES))

LOBJS = $(patsubst %.c, %.lo, $(SOURCES))

LIB_SOURCES = $(SOURCES) $(foreach dir, $(BUILD_DIRS), $(wildcard $(dir)/*.c))

LIB_OBJS = $(patsubst %.c, %.lo, $(LIB_SOURCES))

EXMP_SOURCES = $(wildcard examples/*.c)

TEST_SOURCES = $(wildcard test/*.c)

PROF_SOURCES = $(wildcard profile/*.c)

EXMPS = $(patsubst %.c, %, $(EXMP_SOURCES))

TESTS = $(patsubst %.c, %, $(TEST_SOURCES))

PROFS = $(patsubst %.c, %, $(PROF_SOURCES))

all: $(OBJS) recursive library 

clean:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) clean;)
	rm -f $(OBJS) $(LOBJS) $(TESTS) $(PROFS) $(EXMPS) $(DEFORMATION_LIB) 

distclean: clean
	rm Makefile

profile: all profiler.o
	$(foreach prog, $(PROFS), $(CC) -O2 -std=c99 $(INCS) $(prog).c profiler.o -o $(prog) $(LIBS);)
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) profile;)

recursive:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir);) 

examples: all $(LOBJS) library
	$(foreach prog, $(EXMPS), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(prog) $(LIBS);)

check: $(DEFORMATION_LIB)
ifndef MOD
	$(foreach prog, $(TESTS), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(prog) $(LIBS);)
	$(foreach prog, $(TESTS), $(prog);)
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) check;)
else
	$(foreach dir, $(MOD), $(MAKE) -C $(dir) check;) 
endif

library: library-recursive $(LIB_OBJS)
	$(CC) -fPIC -shared $(LIB_OBJS) $(LIBS2) -o libdeformation.so

library-recursive:
	$(foreach dir, $(BUILD_DIRS), $(MAKE) -C $(dir) library;) 

$(DEFORMATION_LIB): library

install: library
	cp $(DEFORMATION_LIB) $(PREFIX)/lib
	cp *.h $(PREFIX)/include

.PHONY: profile library library-recursive recursive clean check check-recursive all

%.lo: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

%.o: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

BUILD_DIRS = mat_coo mat_csr 

