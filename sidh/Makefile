BUILD_DIR = build
TEST_DIR = test

SOURCES = sidh_elliptic_curve.c\
sidh_isogeny.c\
sidh_private_key.c\
sidh_public_key.c\
sidh_public_param.c\
sidh_quadratic_ext.c\
sidh_shared_key.c\
sidh_util.c\
sidh_public_key_encryption.c\
sidh_elliptic_curve_dlp.c\
sidh_public_key_validation.c

OBJECTS = $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES))

GCC=gcc
GCCFLAGS=-Wall -O3 -march=native -std=c11


.PHONY: clean

all: $(OBJECTS) gen_exec

clean:
	rm -f $(OBJECTS)
	rm -rf $(BUILD_DIR)
	cd $(TEST_DIR) && $(MAKE) clean

$(BUILD_DIR)/%.o:	%.c
	$(GCC) -c $< -o $@ $(GCCFLAGS)

gen_exec: 
	cd $(TEST_DIR) && $(MAKE)

# make the build dir if not exist
$(OBJECTS) : | $(BUILD_DIR)
$(BUILD_DIR) :
	mkdir $(BUILD_DIR)
