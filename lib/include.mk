DIR := lib
DIRS += $(DIR)

_BIO_SRCS := $(shell find $(DIR)/bio -name \*.cc)
_UTIL_SRCS := $(shell find $(DIR)/util -name \*.cc)
_FILESYSTEM_SRCS := $(shell find $(DIR)/filesystem -name \*.cc)
_MATH_SRCS := $(shell find $(DIR)/math -name \*.cc)
_POLYTOPE_SRCS := $(filter-out $(DIR)/polytope/CGALPolytope3D.cc, \
                    $(shell find $(DIR)/polytope -name \*.cc))
_PYTHON_LIBS := $(wildcard $(DIR)/*.py)

_CND_SRCS := $(_BIO_SRCS) $(_UTIL_SRCS) $(_FILESYSTEM_SRCS) $(_MATH_SRCS)

ifdef POLYMAKE_PATH
_CND_SRCS += $(_POLYTOPE_SRCS)
$(_POLYTOPE_SRCS:.cc=$(O)): USE_POLYMAKE = 1
$(_POLYTOPE_SRCS:.cc=$(D)): USE_POLYMAKE = 1
endif

$(DIR)/libcnd$(A): $(_CND_SRCS:.cc=$(O))

SRCS += $(_CND_SRCS)
BINS += 
LIBS += $(DIR)/libcnd$(A)
PYTHON_LIBS += $(_PYTHON_LIBS)

DIST_FILES += $(_BIO_SRCS) $(_UTIL_SRCS) $(_FILESYSTEM_SRCS) \
        $(_POLYTOPE_SRCS) $(_MATH_SRCS) $(DIR)/include.mk \
	$(DIR)/bio/translation/translation.tables $(_PYTHON_LIBS)
