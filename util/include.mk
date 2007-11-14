DIR := util
DIRS += $(DIR)

LOCAL_SRCS := $(filter-out $(DIR)/convexHull2D.cc, $(wildcard $(DIR)/*.cc))
LOCAL_BINS := $(LOCAL_SRCS:.cc=$(E))
LOCAL_PYTHON_SCRIPTS := $(wildcard $(DIR)/*.py)
LOCAL_BASH_SCRIPTS := $(wildcard $(DIR)/*.bash)
SRCS += $(LOCAL_SRCS)
BINS += $(LOCAL_BINS)
PYTHON_SCRIPTS += $(LOCAL_PYTHON_SCRIPTS)
BASH_SCRIPTS += $(LOCAL_BASH_SCRIPTS)
DIST_FILES += $(LOCAL_SRCS) \
              $(LOCAL_PYTHON_SCRIPTS) \
              $(LOCAL_BASH_SCRIPTS) \
              $(DIR)/include.mk
