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

#include "TotalScoreVisitor.hh"
#include "BreakpointSegment.hh"

TotalScoreVisitor::TotalScoreVisitor(Graph& g, Score& totalScore)
	: vertexColors(get(vertex_color, g)),
	  vertexIndices(get(vertex_index, g)),
	  breakpointSegments(get(breakpoint_segment_t(), g)),
	  alignedSegments(get(aligned_segments_t(), g)),
	  totalScore(totalScore) {
}

void TotalScoreVisitor::start_vertex(Vertex u, const Graph& g) {
	BreakpointSegment* seg = breakpointSegments[u];
	Score maxScore = std::numeric_limits<Score>::min();
	for (size_t index = 0; index < seg->getNumIndices(); ++index) {
		if (seg->getScore(index) > maxScore) {
			maxScore = seg->getScore(index);
		}
	}
	totalScore += maxScore;
}
