DIR := include
_BIO_SRCS := $(shell find $(DIR)/bio -name \*.hh)
_POLYTOPE_SRCS := $(shell find $(DIR)/polytope -name \*.hh)
_UTIL_SRCS := $(shell find $(DIR)/util -name \*.hh)
_MATH_SRCS := $(shell find $(DIR)/math -name \*.hh)
_FILESYSTEM_SRCS := $(shell find $(DIR)/filesystem -name \*.hh) \
                    $(DIR)/filesystem.hh
DIST_FILES += $(_BIO_SRCS) $(_UTIL_SRCS) $(_POLYTOPE_SRCS) $(_MATH_SRCS) \
              $(_FILESYSTEM_SRCS) $(DIR)/include.mk
DIST_DIRS += $(DIR)/boost
