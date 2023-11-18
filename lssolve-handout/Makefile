CC=gcc
CPPFLAGS=-Imsptools/include
CFLAGS=-Wall -Werror=vla -std=c11
LDFLAGS=-Lmsptools/lib
LDLIBS=-lmsptools -lm

OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
	LDLIBS=-lopenblas -lmsptools -lm
endif
ifeq ($(OS_NAME),Darwin)
	LDLIBS=-llapack -lblas -lmsptools -lm
endif

.PHONY: all clean test1 test2 test

all: lssolve test_call_dgels lssolve-handin.tar

lssolve: call_dgels.o msptools/lib/libmsptools.a

test_call_dgels: call_dgels.o msptools/lib/libmsptools.a

test: test1 test2

test1: test_call_dgels
	@echo "****** test1 ******"
	@./test_call_dgels 

test2: lssolve
	@echo "****** test2 ******"
	@./test_lssolve.sh 

lssolve-handin.tar: call_dgels.c lssolve.c
	tar cvf lssolve-handin.tar $^

msptools/lib/libmsptools.a:
	$(MAKE) -C msptools

clean:
	-$(RM) -rf lssolve test_call_dgels *.o lssolve-handin.tar data/x*.txt
