DIR := apps/mercator
LOCAL_SRCS := $(wildcard $(DIR)/*.cc)
LOCAL_HEADERS := $(wildcard $(DIR)/*.hh)
LOCAL_BINS := $(foreach bin,mercator,$(DIR)/$(bin)$(E))
SRCS += $(LOCAL_SRCS)
BINS += $(LOCAL_BINS)
DIST_FILES += $(LOCAL_SRCS) $(LOCAL_HEADERS) $(DIR)/include.mk $(DIR)/README

MERCATOR_BUILD_VERSION := 0.4
BUILD_DATE := $(shell date +"%F %R")

$(LOCAL_BINS): $(LOCAL_SRCS:.cc=$(O))
$(LOCAL_BINS): CPPFLAGS += -D BUILD_VERSION="\"$(MERCATOR_BUILD_VERSION)\"" -D BUILD_DATE="\"$(BUILD_DATE)\""

include $(wildcard $(DIR)/*/include.mk)
