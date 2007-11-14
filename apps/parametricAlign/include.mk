DIR := apps/parametricAlign
LOCAL_SRCS := $(wildcard $(DIR)/*.cc)
LOCAL_HEADERS := $(wildcard $(DIR)/*.hh)
LOCAL_MAINS := faAlignPolytope lengthEffect vertexAlignments vertexParams faAlignPolytopeRow faAlignSummary faAlignSummarySpeedTest normalFan facetNormals points2polytope minkowskiSum processAlignmentPolytope findMaxIdentityAlignments maxIdentityAlign checkPointsInFacets summarizeAlignment
LOCAL_BINS := $(foreach bin,$(LOCAL_MAINS),$(DIR)/$(bin)$(E))

DIST_FILES += $(LOCAL_SRCS) $(LOCAL_HEADERS) $(DIR)/include.mk
$(LOCAL_BINS): $(DIR)/polyalign$(O) \
               $(DIR)/VectorMultiSet$(O) \
               $(DIR)/AnnotatedPolytope$(O)
$(LOCAL_BINS): USE_POLYMAKE = 1
$(LOCAL_BINS): USE_QHULL = 1

ifdef POLYMAKE_PATH
SRCS += $(LOCAL_SRCS)
BINS += $(LOCAL_BINS)
$(LOCAL_SRCS:.cc=$(D)): USE_POLYMAKE = 1
$(LOCAL_SRCS:.cc=$(D)): USE_QHULL = 1
else
$(warning parametricAlign will not be built because POLYMAKE_PATH environment variable is not set)
endif
