###############################################################################
## Copyright (C) 2010 Sebastian Pancratz                                     ##
###############################################################################

CC = gcc
CFLAGS = -g -ansi -Wall -O2
LIBS = -lgmp -lm
INCS = -I ../ -I ../monv/ -I ../basis/ -I ../padics/ -I ../padics/padic/

OBJS = $(patsubst %.c, %.o, $(wildcard *.c))

all: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c $(INCS) $< -o $@

clean: 
	rm *.o

