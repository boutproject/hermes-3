
BOUT_TOP	= ../..

TARGET = hermes-3

DIRS = src

SOURCEC		= hermes-3.cxx

# Capture the git version, to be printed in the outputs
GIT_VERSION := $(shell git describe --abbrev=40 --dirty --always --tags)
CXXFLAGS += -DHERMES_VERSION=\"$(GIT_VERSION)\"

include $(BOUT_TOP)/make.config

check-unit-tests: src.a
	$(MAKE) -C tests/unit check

check-integrated-tests: hermes-3
	@cd tests/integrated; ./test_suite

check: check-unit-tests check-integrated-tests

clean::
	make -C tests/unit clean
