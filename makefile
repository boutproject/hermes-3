
BOUT_TOP	= ../..

TARGET = hermes-2

DIRS = atomicpp

SOURCEC		= hermes-2.cxx div_ops.cxx loadmetric.cxx radiation.cxx neutral-model.cxx diffusion2d.cxx recycling.cxx full-velocity.cxx mixed.cxx

include $(BOUT_TOP)/make.config
