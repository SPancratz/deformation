SOURCES = $(wildcard *.c)

OBJS = $(patsubst %.c, %.o, $(SOURCES))

LIB_OBJS = $(patsubst %.c, %.lo, $(SOURCES))

TEST_SOURCES = $(wildcard test/*.c)

PROF_SOURCES = $(wildcard profile/*.c)

TESTS = $(patsubst %.c, %, $(TEST_SOURCES))

PROFS = $(patsubst %.c, %, $(PROF_SOURCES))

all: $(OBJS)

library: $(LIB_OBJS)

profile:
	$(foreach prog, $(PROFS), $(CC) -O2 -std=c99 $(INCS) $(prog).c ../profiler.o -o $(prog) $(LIBS) -ldeformation;)
        
%.o: %.c
	$(CC) $(CFLAGS) -c $(INCS) $< -o $@

%.lo: %.c
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -f $(OBJS) $(LIB_OBJS) $(TESTS)	

check: library
	$(foreach prog, $(TESTS), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(prog) $(LIBS) -ldeformation;)
	$(foreach prog, $(TESTS), valgrind --leak-check=full $(prog);)

.PHONY: profile clean check all
