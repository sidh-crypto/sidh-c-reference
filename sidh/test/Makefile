INCLUDE_DIRS = -I.. -I/usr/local/include
LIB_DIRS = -L/usr/local/lib
BUILD_DIR = ../build
OBJECTS = $(wildcard $(BUILD_DIR)/*)
SOURCES = $(wildcard test*.c)
TESTS = $(patsubst %.c, %, $(SOURCES))
EXECS = $(patsubst %.c, %, $(SOURCES))

GCC = gcc
GCCFLAGS = -Wall -O3 -march=native -std=c11
LIBS = $(LIB_DIRS) -lgmp -lm

.PHONY: clean

all:	$(TESTS)

clean:
	rm -f $(EXECS)

test_%:	$(OBJECTS)
	if test -e test_$*.c; then\
		$(GCC) -o test_$* test_$*.c $(OBJECTS) $(GCCFLAGS) $(INCLUDE_DIRS) $(LIBS);\
	fi
