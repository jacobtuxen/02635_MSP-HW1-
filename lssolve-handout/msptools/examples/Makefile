CC=gcc
CPPFLAGS=-I../include
CFLAGS=-Wall -g -std=c99
LDFLAGS=-L../lib
LDLIBS=-lm -lmsptools
EXAMPLES=$(basename $(wildcard *example*.c))

.PHONY: all clean

all: $(EXAMPLES)

$(EXAMPLES):  ../lib/libmsptools.a

clean:
	-$(RM) $(EXAMPLES)	
	-$(RM) -r *example*.dSYM
