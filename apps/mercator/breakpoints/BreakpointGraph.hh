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

#ifndef __BREAKPOINT_GRAPH_HH__
#define __BREAKPOINT_GRAPH_HH__

#include "util/io.hh"

#include <iostream>
#include <vector>
#include <stdexcept>
#include <utility>

// General utilities
#include "boost/bind.hpp"
#include "boost/property_map/property_map.hpp"
// Graph utilities
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/kruskal_min_spanning_tree.hpp"
#include "boost/graph/graph_utility.hpp"
#include "boost/graph/undirected_dfs.hpp"
using namespace boost;

#include "Genome.hh"
#include "Segment.hh"
#include "SegmentSet.hh"
#include "BreakpointSegment.hh"
#include "AlignedSegments.hh"
#include "AlignedSegmentsWeight.hh"
#include "TotalScoreVisitor.hh"
#include "ScoringVisitor.hh"
#include "TracebackVisitor.hh"
#include "graph_types.hh"

#include "bio/sdb.hh"
#include "util/io.hh"
#include "bio/alignment/DNAScoringMatrix.hh"
#include "bio/alignment/AmbiguousDNAScoringMatrix.hh"
#include "util/matrix.hh"
#include "util/string.hh"
#include "util/options.hh"
using bio::alignment::DNAScoringMatrix;
using bio::alignment::AmbiguousDNAScoringMatrix;

class BreakpointGraph {
public:
	BreakpointGraph();
	void readHomologyMap(std::istream& strm);

	void addAllEdges(util::Matrix<double>& distances,
					 bool removeColinearEdges = false);
	void addSomeEdges(std::istream& strm, util::Matrix<double>& distances);

	void calcMidpointScore(util::Matrix<double>& distances);
	
	void calcBestBreakpoints(util::Matrix<double>& distances,
							 size_t maxSegmentPositions);

	void setBreakpoints(std::istream& stream);
	void refineMap();
	
	void writeBreakpointSegments(std::ostream& strm);
	void writeMSTEdges(std::ostream& strm);
	void writeBreakpoints(std::ostream& strm);
	void writeMap(std::ostream& strm);

	void writeAlignedSegments(std::ostream& stream);
	
private:
	std::vector<Segment*> segments;
	std::vector<SegmentSet*> segmentSets;
	std::vector<BreakpointSegment*> breakpointSegments;
	std::vector<AlignedSegments*> alignments;
	Graph g;
	
	void makeGraph();
	void loadSegmentSets(std::istream& strm);
	void loadSegmentSet(const std::string& line);
	void sortSegments();
	Segment* makeStartSegment(Genome* genome,
							  const std::string& chrom) const;
	Segment* makeEndSegment(Genome* genome,
							const std::string& chrom) const;
	void addChromosomeEndSegments();
	void makeBreakpointSegments();
	void makeAlignedSegments();
	void removeColinearAlignedSegments();
	
	void readAlignments(util::Matrix<double>& distances);
	void resetScores();
	void initBounds();
	void updateBounds(size_t minUpdateLength);
	void setBestBreakpoints();
	size_t maxSegmentBoundedLength();
	void updateSegmentScores();
};

#endif // __BREAKPOINT_GRAPH_HH__
