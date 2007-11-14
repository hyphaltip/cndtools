DIR := apps
DIST_FILES += $(DIR)/include.mk
include $(wildcard $(DIR)/*/include.mk)
