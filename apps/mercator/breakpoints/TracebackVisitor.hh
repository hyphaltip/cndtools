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

#ifndef __TRACEBACK_VISITOR_HH__
#define __TRACEBACK_VISITOR_HH__

#include "graph_types.hh"

class TracebackVisitor : public default_dfs_visitor {
private:
	VertexColorMap vertexColors;
	VertexIndexMap vertexIndices;
	BreakpointSegmentMap breakpointSegments;
	AlignedSegmentsMap alignedSegments;

	void traceback(BreakpointSegment* seg1,
				   BreakpointSegment* seg2,
				   AlignedSegments* alignment);
	
public:	
	TracebackVisitor(Graph& g);
	
	void start_vertex(Vertex u, const Graph& g);
	
	void discover_vertex(Vertex u, const Graph& g);
};

#endif // __TRACEBACK_VISITOR_HH__
