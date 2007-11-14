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

#include "ScoringVisitor.hh"
#include "BreakpointSegment.hh"
#include "AlignedSegments.hh"

ScoringVisitor::ScoringVisitor(Graph& g)
	: vertexColors(get(vertex_color, g)),
	  vertexIndices(get(vertex_index, g)),
	  breakpointSegments(get(breakpoint_segment_t(), g)),
	  alignedSegments(get(aligned_segments_t(), g)) {
}
	
void ScoringVisitor::addScores(BreakpointSegment* seg1,
							   BreakpointSegment* seg2,
							   AlignedSegments* alignment) {
	bool flip = (alignment->getSegment1() == seg2);
		
	for (size_t index1 = 0; index1 < seg1->getNumIndices(); ++index1) {
		Score maxScore = std::numeric_limits<Score>::min();
		for (size_t index2 = 0; index2 < seg2->getNumIndices(); ++index2) {
			Score score = seg2->getScore(index2);
			score += (flip ?
					  alignment->score(index2, index1) :
					  alignment->score(index1, index2));
			if (score > maxScore) {
				maxScore = score;
			}
		}
		seg1->setScore(index1, maxScore + seg1->getScore(index1));
	}
}
	
void ScoringVisitor::finish_vertex(Vertex u, const Graph& g) {
	//std::cerr << "Finishing vertex " << vertexIndices[u] << '\n';
		
	OutEdgeIterator out_i, out_end;
	for (tie(out_i, out_end) = out_edges(u, g); out_i != out_end; ++out_i) {
		Vertex v = target(*out_i, g);
		if (vertexColors[v] == Color::black()) {
			addScores(breakpointSegments[u],
					  breakpointSegments[v],
					  alignedSegments[*out_i]);
		}
	}
}
