SHELL=/bin/sh

.SUFFIXES:
.SUFFIXES: .cc .o .d .hh .a

.LIBPATTERNS = lib%.a lib%.so

DIST_NAME = cndsrc
DIST_DATE = $(shell date +"%Y.%m.%d")
DIST_PREFIX = $(DIST_NAME)-$(DIST_DATE)

prefix ?= /usr/local
exec_prefix ?= $(prefix)
bindir ?= $(exec_prefix)/bin
includedir ?= $(prefix)/include
libdir ?= $(exec_prefix)/lib

ifeq ($(OS),MACOSX)
  CXXOPT = -O3 -ffast-math #-fast -mcpu=G4
else
  CXXOPT = -O3 -ffast-math
endif

CXXDEBUG = -ggdb #-D_GLIBCXX_DEBUG

CC = $(CXX)
CXXFLAGS += -Wall
CPPFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE 
CPPFLAGS += -Iinclude
LDFLAGS += #-ggdb

ifeq ($(OS),MACOSX)
  CXXFLAGS += -Wno-long-double
endif

vpath lib%.a lib

ifeq ($(DEBUG),)
  CPPFLAGS += -DNDEBUG
  CXXFLAGS += $(CXXOPT) # -ggdb
else
  ifeq ($(DEBUG),profile)
    SUFFIX :=-p
    CXXFLAGS += -pg $(CXXOPT)
    LDFLAGS += -pg
    ifeq ($(PROFILE),detailed)
      CXXFLAGS += -fno-inline -fno-default-inline
      LDFLAGS += -Wl,-z,muldefs
    endif
  else
    SUFFIX :=-d
    CXXFLAGS += $(CXXDEBUG)
    LDFLAGS += $(CXXDEBUG)

    ifeq ($(DEBUG),opt)
      SUFFIX :=$(SUFFIX)$(patsubst -%,%,$(CXXOPT))
    endif
  endif
endif

E := $(SUFFIX)
O := $(SUFFIX).o
A := $(SUFFIX).a
D := .d

# CGAL SUPPORT
ifneq ($(CGAL_MAKEFILE),)
  include $(CGAL_MAKEFILE)
  CXXFLAGS += $(CGAL_CXXFLAGS)
  LDLIBS += $(CGAL_LDFLAGS)
endif

# POLYMAKE SUPPORT
# The environment variable POLYMAKE_PATH must be set appropriately
CPPFLAGS += $(if $(USE_POLYMAKE),$(addprefix -I$(POLYMAKE_PATH)/, \
                                             apps/polytope/include \
                                             lib/poly_client/include \
                                             lib/PTL/include \
                                             lib/PTL/include/std \
                                             lib/gmp_wrapper/include))
CPPFLAGS += $(if $(USE_POLYMAKE),-DALL_PLAUSIBLE_CHECKS=$(if $(DEBUG),1 -DPOLY_DEBUG,0))
CXXFLAGS += $(if $(USE_POLYMAKE),-ftemplate-depth-200)
LDLIBS += $(if $(USE_POLYMAKE),$(foreach lib,polytope poly,-l$(lib)) -lgmp)
LDFLAGS += $(if $(USE_POLYMAKE),-L$(POLYMAKE_PATH)/lib)

# Support for Christophe Weibel's MINKSUM software
CPPFLAGS += $(if $(USE_MINKSUM),-I$(MINKSUM_PATH)/include)
LDLIBS += $(if $(USE_MINKSUM),$(foreach lib,MINKSUM cdd gmp wrapgmp-gcc3, -l$(lib)))
LDFLAGS += $(if $(USE_MINKSUM),-L$(MINKSUM_PATH)/lib)

# Support for Qhull
CPPFLAGS += $(if $(USE_QHULL),-I$(QHULL_PATH)/include)
LDLIBS += $(if $(USE_QHULL), -lqhull)
LDFLAGS += $(if $(USE_QHULL),-L$(QHULL_PATH)/lib)

# Cancel implicit rule for making binary directly from source.  Prefer
# to build all objects first so that dependencies are followed
# correctly.
%: %.cc

# Dependency file creation rule
%$(D): %.cc
	@$(CXX) -MM -MT '$(*)$$(O) $@' $(CPPFLAGS) $< > $@

%$(O): %.cc
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

%$(A):
	rm -f $@
	$(AR) qc $@ $^

DIRS :=
SRCS :=
LIBS :=
BINS :=
DIST_FILES := makefile README COPYING
DIST_DIRS :=
PYTHON_LIBS :=
PYTHON_SCRIPTS :=
BASH_SCRIPTS :=

OBJS = $(SRCS:.cc=$(O))
DEPENDS = $(SRCS:.cc=$(D))

# Make sure that "all" is the default target before including other makefiles
all:

# Recursively include submakefiles
include $(wildcard */include.mk)

# Include dependency files for all targets except for clean
ifeq (,$(findstring clean,$(MAKECMDGOALS)))
-include $(DEPENDS)
endif

.PHONY: all clean install depends $(DIRS)

# Make all binaries depend on the library
$(BINS): -lcnd$(SUFFIX)

ALL_TARGETS = $(LIBS) $(BINS)

all: $(ALL_TARGETS)

depends: $(DEPENDS)

$(bindir):
	mkdir -p $(bindir)

install: $(bindir) $(BINS)
	install -sp $(addsuffix $(EXE),$(BINS)) $(bindir)
	for p in $(PYTHON_SCRIPTS); do cp $$p $(bindir)/`basename $$p .py`; chmod +x $(bindir)/`basename $$p .py`; done
	for p in $(PYTHON_LIBS); do cp $$p $(bindir); done
	for b in $(BASH_SCRIPTS); do cp $$b $(bindir)/`basename $$b .bash`; chmod +x $(bindir)/`basename $$b .bash`; done

dist:
	rm -rf $(DIST_PREFIX)
	mkdir $(DIST_PREFIX)
	tar hczf $(DIST_PREFIX)/$(DIST_PREFIX).tar.gz $(DIST_FILES) $(DIST_DIRS)
	cd $(DIST_PREFIX); tar zxf $(DIST_PREFIX).tar.gz; rm $(DIST_PREFIX).tar.gz
	tar czf $(DIST_PREFIX).tar.gz $(DIST_PREFIX)
	rm -rf $(DIST_PREFIX)

clean:
	rm -f $(DEPENDS) $(OBJS) $(BINS) $(LIBS)

cleandepends:
	rm -f $(DEPENDS)

$(DIRS):
	$(MAKE) $(filter $@/%, $(ALL_TARGETS))
