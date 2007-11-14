DIR := apps/mercator/util
LOCAL_SRCS := $(wildcard $(DIR)/*.cc)
LOCAL_HEADERS := $(wildcard $(DIR)/*.hh)

LOCAL_MAINS := anchors2fa gff2anchors makeAlignmentInput matchSegment \
	       maskRepetitive
LOCAL_OBJS := $(LOCAL_SRCS:.cc=$(O))
LOCAL_MAIN_OBJS := $(foreach bin,$(LOCAL_MAINS),$(DIR)/$(bin)$(O))
LOCAL_SHARED_OBJS := $(filter-out $(LOCAL_MAIN_OBJS), $(LOCAL_OBJS))
LOCAL_BINS := $(foreach bin,$(LOCAL_MAINS),$(DIR)/$(bin)$(E))

$(LOCAL_BINS): $(LOCAL_SHARED_OBJS)

LOCAL_PYTHON_SCRIPTS := $(wildcard $(DIR)/*.py)
LOCAL_BASH_SCRIPTS := $(wildcard $(DIR)/*.bash)

PYTHON_SCRIPTS += $(LOCAL_PYTHON_SCRIPTS)
BASH_SCRIPTS += $(LOCAL_BASH_SCRIPTS)

SRCS += $(LOCAL_SRCS)
BINS += $(LOCAL_BINS)
DIST_FILES += $(LOCAL_SRCS) $(LOCAL_HEADERS) $(DIR)/include.mk \
              $(LOCAL_PYTHON_SCRIPTS) $(LOCAL_BASH_SCRIPTS)
