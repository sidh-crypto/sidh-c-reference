BUILD_DIR = ../build
OBJECTS = $(wildcard $(BUILD_DIR)/*)
SOURCES = $(wildcard test*.c)
TESTS = $(patsubst %.c, %, $(SOURCES))
INCLUDE_DIRS = -I.. -I/usr/local/include
LIB_DIRS = -L$(BUILD_DIR) -L/usr/local/lib
EXECS = $(patsubst %.c, %, $(SOURCES))

GCC = gcc
GCCFLAGS = -Wall -O0 -g -march=native -std=c11 -static
LIBS = $(LIB_DIRS) -lsidh -lgmp -lm

.PHONY: clean

all:	$(TESTS)

clean:
	rm -f $(EXECS)

test_%:
	$(GCC) -o test_$* test_$*.c $(GCCFLAGS) $(INCLUDE_DIRS) $(LIBS)
