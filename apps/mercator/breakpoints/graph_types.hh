/* Copyright (c) 2006
   Colin Dewey (University of Wisconsin-Madison)
   cdewey@biostat.wisc.edu
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __GRAPH_TYPES_HH__
#define __GRAPH_TYPES_HH__

#include "boost/property_map/property_map.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/undirected_dfs.hpp"
using namespace boost;

#include "types.hh"
#include "AlignedSegmentsWeight.hh"

struct breakpoint_segment_t {
	typedef vertex_property_tag kind;
};

struct aligned_segments_t { 
	typedef edge_property_tag kind;
};

typedef property<breakpoint_segment_t, BreakpointSegment*,
				 property<vertex_color_t, int> > VertexProperty;
typedef property<edge_index_t, size_t,
				 property<edge_color_t, int,
						  property<edge_weight_t, AlignedSegmentsWeight,
								   property<aligned_segments_t, AlignedSegments*> > > > EdgeProperty;

typedef adjacency_list<vecS,
					   vecS,
					   undirectedS,
					   VertexProperty,
					   EdgeProperty> Graph;

typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::out_edge_iterator OutEdgeIterator;

typedef property_map<Graph, edge_index_t>::type EdgeIndexMap;
typedef property_map<Graph, vertex_index_t>::type VertexIndexMap;
typedef property_map<Graph, vertex_color_t>::type VertexColorMap;
typedef property_map<Graph, edge_color_t>::type EdgeColorMap;
typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap;
typedef property_map<Graph, breakpoint_segment_t>::type BreakpointSegmentMap;
typedef property_map<Graph, aligned_segments_t>::type AlignedSegmentsMap;

typedef property_traits<VertexColorMap>::value_type ColorValue;
typedef color_traits<ColorValue> Color;

#endif // __GRAPH_TYPES_HH__
