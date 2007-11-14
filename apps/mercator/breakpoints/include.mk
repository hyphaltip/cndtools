DIR := apps/mercator/breakpoints
LOCAL_MAINS := makeBreakpointGraph findBreakpoints breakMap \
               makeBreakpointAlignmentInput breakpointMultipleAlignments
LOCAL_MAIN_SRCS := $(foreach bin,$(LOCAL_MAINS),$(DIR)/$(bin).cc)
LOCAL_SRCS := $(wildcard $(DIR)/*.cc)
LOCAL_HEADERS := $(wildcard $(DIR)/*.hh)
LOCAL_BINS := $(foreach bin,$(LOCAL_MAINS),$(DIR)/$(bin)$(E))
LOCAL_OBJS := $(patsubst %.cc,%$(O),$(filter-out $(LOCAL_MAIN_SRCS),$(LOCAL_SRCS)))

SRCS += $(LOCAL_SRCS)
BINS += $(LOCAL_BINS)
DIST_FILES += $(LOCAL_SRCS) $(LOCAL_HEADERS) $(DIR)/include.mk

$(LOCAL_BINS): $(LOCAL_OBJS)

#include $(wildcard $(DIR)/*/include.mk)
