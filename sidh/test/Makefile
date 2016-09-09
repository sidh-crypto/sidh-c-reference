INC_DIR = ..
BUILD_DIR = ../build
OBJS = $(wildcard $(BUILD_DIR)/*)
SOURCES = $(wildcard test*.c)
TESTS = $(patsubst %.c, %, $(SOURCES))
EXEC = $(patsubst %.c, %, $(SOURCES))

GCC = gcc
GCCFLAGS = -Wall -O3 -march=native -std=c11 -I$(INC_DIR)
LIBS = -lgmp -lm

.PHONY: clean

all:	$(TESTS)

clean:
	rm -f $(EXEC)

test_%:	$(OBJS)
	if test -e test_$*.c; then\
		$(GCC) -o test_$* test_$*.c $(OBJS) $(GCCFLAGS) $(LIBS);\
	fi
