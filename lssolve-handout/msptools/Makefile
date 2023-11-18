CC=gcc
CPPFLAGS=-Iinclude
CFLAGS=-Wall -O3 -std=c99
LDFLAGS=
LDLIBS=-lm

INSTALL_PREFIX=/usr/local
objs=$(patsubst %.c,%.o,$(wildcard src/*.c))

.PHONY: all examples test clean

all: examples 

src/%.o: src/%.c

lib/libmsptools.a: $(objs)
	mkdir -p lib
	ar rcs $@ $^

examples: lib/libmsptools.a
	$(MAKE) --directory=examples

test: lib/libmsptools.a
	$(MAKE) --directory=tests

install: lib/libmsptools.a
	mkdir -p $(INSTALL_PREFIX)/lib
	mkdir -p $(INSTALL_PREFIX)/include/msptools
	cp lib/libmsptools.a $(INSTALL_PREFIX)/lib
	cp include/*.h $(INSTALL_PREFIX)/include/msptools

clean:
	-$(RM) src/*.o lib/libmsptools.a
	-$(RM) data/*_copy.txt
	$(MAKE) --directory=examples clean
	$(MAKE) --directory=tests clean

