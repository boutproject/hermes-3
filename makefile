
BOUT_TOP	= ../..

TARGET = hermes-3

DIRS = src

SOURCEC		= hermes-3.cxx

include $(BOUT_TOP)/make.config

check-unit-tests:
	$(MAKE) -C unit_tests check

check: check-unit-tests
