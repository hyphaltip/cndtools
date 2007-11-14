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

#include "BreakpointGraph.hh"
#include "util/io/line/InputStream.hh"

BreakpointGraph::BreakpointGraph()
	: segments(),
	  segmentSets(),
	  breakpointSegments(),
	  alignments(),
	  g()
{
}

void BreakpointGraph::readHomologyMap(std::istream& strm) {
		std::cerr << "Reading segments from homology map...";
		loadSegmentSets(strm);
		std::cerr << "done.\n";
		std::cerr << "Homology map: "
				  << segments.size() << " segments in "
				  << segmentSets.size() << " homology sets" << '\n';

		std::cerr << "Sorting segments...";
		sortSegments();
		std::cerr << "done.\n";

		std::cerr << "Adding zero-length chromosome end segments...";
		addChromosomeEndSegments();
		std::cerr << "done.\n";
		
		std::cerr << "Making breakpoint segments...";
		makeBreakpointSegments();
		std::cerr << "done.\n";

		std::cerr << "Creating breakpoint graph...";
		makeGraph();
		std::cerr << "done.\n";
}

void BreakpointGraph::makeGraph() {
	g = Graph(breakpointSegments.size());
	// Add breakpoint segments as property of vertices
	BreakpointSegmentMap breakpoint_segments = get(breakpoint_segment_t(), g);
	for (size_t i = 0; i < breakpointSegments.size(); ++i) {
		breakpoint_segments[i] = breakpointSegments[i];
	}
}

void BreakpointGraph::removeColinearAlignedSegments() {
	alignments.erase(std::remove_if(alignments.begin(),
									alignments.end(),
									std::mem_fun(&AlignedSegments::isColinear)),
					 alignments.end());
}

void BreakpointGraph::addAllEdges(util::Matrix<double>& distances,
								  bool removeColinearEdges) {
	std::cerr << "Connecting aligned segments...";
	makeAlignedSegments();
	std::cerr << "done.\n";

	if (removeColinearEdges) {
		std::cerr << "Removing colinear aligned segments...";
		removeColinearAlignedSegments();
		std::cerr << "done.\n";
	}
	
	for (std::vector<AlignedSegments*>::const_iterator it = alignments.begin();
		 it != alignments.end(); ++it) {
		AlignedSegments* as = *it;
		Edge e;
		bool success;
		tie(e, success) = add_edge((*it)->getSegment1()->getNum(),
								   (*it)->getSegment2()->getNum(),
								   g);

		if (not success) {
			std::cerr << "Could not add edge\n";
		}

		size_t genomeNum1 = as->getSegment1()->getGenome()->getNum();
		size_t genomeNum2 = as->getSegment2()->getGenome()->getNum();
		double distance = distances(genomeNum1, genomeNum2);

		// Initialize properties for this edge
		put(aligned_segments_t(), g, e, as);
		put(edge_weight_t(), g, e, AlignedSegmentsWeight(as, distance));
	}
}

void BreakpointGraph::addSomeEdges(std::istream& stream,
								   util::Matrix<double>& distances) {
	size_t edgeNum;
	size_t seg1Num;
	size_t seg2Num;
	std::string strand1;
	std::string strand2;
	
	while (stream >> edgeNum
		   >> seg1Num >> seg2Num >> strand1 >> strand2) {
		AlignedSegments* as = new AlignedSegments(breakpointSegments[seg1Num],
												  breakpointSegments[seg2Num],
												  strand1 == "+",
												  strand2 == "+");
		alignments.push_back(as);
		Edge e;
		bool success;
		tie(e, success) = add_edge(as->getSegment1()->getNum(),
								   as->getSegment2()->getNum(),
								   edgeNum,
								   g);

		size_t genomeNum1 = as->getSegment1()->getGenome()->getNum();
		size_t genomeNum2 = as->getSegment2()->getGenome()->getNum();
		double distance = distances(genomeNum1, genomeNum2);

		// Initialize properties for this edge
		put(aligned_segments_t(), g, e, as);
		put(edge_weight_t(), g, e, AlignedSegmentsWeight(as, distance));
	}
}

void BreakpointGraph::setBreakpoints(std::istream& stream) {
	size_t segNum;
	size_t breakpoint;
	
	while (stream >> segNum >> breakpoint) {
		breakpointSegments[segNum]->setRelBreakpoint(breakpoint);
	}
}

void BreakpointGraph::loadSegmentSets(std::istream& strm) {
	util::io::line::InputStream lineStream(strm);
	std::string line;
	while (lineStream >> line) {
		loadSegmentSet(line);
	}
}

void BreakpointGraph::loadSegmentSet(const std::string& line) {
	util::string::Converter<size_t> StringToNum;
	std::vector<std::string> tokens;
	util::string::split(line, std::back_inserter(tokens), "\t");
	if (tokens.size() % 5 != 1) {
		throw std::runtime_error("Invalid homology map line: " + line);
	}
	SegmentSet* set = new SegmentSet();
	set->setNum(StringToNum(tokens[0]));
	for (size_t i = 1; i < tokens.size(); i += 5) {
		Segment* seg = new Segment(Genome::getGenome(tokens[i], true),
								   tokens[i + 1],
								   StringToNum(tokens[i + 2]),
								   StringToNum(tokens[i + 3]),
								   tokens[i + 4].at(0));
		set->addSegment(seg);
		segments.push_back(seg);
	}
	segmentSets.push_back(set);
}	

struct SegmentSorter {
	bool operator()(const Segment* x, const Segment* y) const {
		return *x < *y;
	}
};

void BreakpointGraph::sortSegments() {
	std::sort(segments.begin(), segments.end(), SegmentSorter());
}

Segment* BreakpointGraph::makeStartSegment(Genome* genome,
										   const std::string& chrom) const {
	return new Segment(genome, chrom, 0, 0, '+');
}

Segment* BreakpointGraph::makeEndSegment(Genome* genome,
										 const std::string& chrom) const {
	size_t chromLen = genome->getChromLen(chrom);
	return new Segment(genome, chrom, chromLen, chromLen, '+');
}

// Precondition: segments must be sorted
// Postcondition: segments remain sorted
void BreakpointGraph::addChromosomeEndSegments() {
	std::vector<Segment*> newSegments;
	Segment* prevSeg = NULL;
	std::vector<Segment*>::const_iterator pos;
	for (pos = segments.begin(); pos != segments.end(); ++pos) {
		Segment* currSeg = *pos;
		if (prevSeg == NULL or not currSeg->isOnSameChrom(*prevSeg)) {
			if (prevSeg != NULL) {
				newSegments.push_back(makeEndSegment(prevSeg->getGenome(),
													 prevSeg->getChrom()));
			}
			newSegments.push_back(makeStartSegment(currSeg->getGenome(),
												   currSeg->getChrom()));
		}
		newSegments.push_back(currSeg);
		prevSeg = currSeg;
	}

	if (prevSeg != NULL) {
		newSegments.push_back(makeEndSegment(prevSeg->getGenome(),
											 prevSeg->getChrom()));
	}

	segments.swap(newSegments);
}

// Precondition: segments must be sorted
void BreakpointGraph::makeBreakpointSegments() {
	Segment* prevSeg = NULL;
	std::vector<Segment*>::const_iterator pos;
	for (pos = segments.begin(); pos != segments.end(); ++pos) {
		Segment* currSeg = *pos;
		if (prevSeg != NULL && prevSeg->isOnSameChrom(*currSeg)) {
			breakpointSegments.push_back(new BreakpointSegment(prevSeg, currSeg));
// 			std::cerr << "Made breakpoint segment: "
// 					  << *breakpointSegments.back() << '\n'
// 					  << "From " << *prevSeg << " and " << *currSeg << '\n';
		}
		prevSeg = currSeg;
	}

	// Number segments
	for (size_t num = 0; num < breakpointSegments.size(); ++num) {
		breakpointSegments[num]->setNum(num);
	}
}

struct AlignedSegmentsComparer {
	bool operator()(const AlignedSegments* as1,
					const AlignedSegments* as2) const {
		return (as1->getSegment1() < as2->getSegment1())
			or (as1->getSegment1() == as2->getSegment1()
				and as1->getSegment2() < as2->getSegment2());
	}
};

void BreakpointGraph::makeAlignedSegments() {
	std::for_each(segmentSets.begin(), segmentSets.end(),
				  boost::bind(&SegmentSet::makeAlignedSegments, _1,
							  ref(alignments)));
}

void BreakpointGraph::writeMSTEdges(std::ostream& stream) {
	std::vector<Edge> edges;
	kruskal_minimum_spanning_tree(g, std::back_inserter(edges));

	BreakpointSegmentMap segmentsMap = get(breakpoint_segment_t(), g);
	AlignedSegmentsMap alignedSegmentsMap = get(aligned_segments_t(), g);	
	for (size_t i = 0; i < edges.size(); ++i) {
		AlignedSegments& as = *alignedSegmentsMap[edges[i]];
		stream << i << '\t'
			   << segmentsMap[source(edges[i], g)]->getNum() << '\t'
			   << segmentsMap[target(edges[i], g)]->getNum() << '\t'
			   << (as.isPrefix1() ? '+' : '-') << '\t'
			   << (as.isPrefix2() ? '+' : '-') << '\n';
	}
}

void BreakpointGraph::readAlignments(util::Matrix<double>& distances) {
	AlignedSegmentsMap alignedSegmentsMap = get(aligned_segments_t(), g);
	EdgeIndexMap edgeIndexMap = get(edge_index_t(), g);	
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		alignedSegmentsMap[*ei]->readAlignment(edgeIndexMap[*ei], distances);
	}
}

void BreakpointGraph::initBounds() {
	std::for_each(breakpointSegments.begin(), breakpointSegments.end(),
				  std::mem_fun(&BreakpointSegment::initBounds));
}

void BreakpointGraph::updateBounds(size_t minUpdateLength) {
	for (size_t i = 0; i < breakpointSegments.size(); ++i) {
// 		if (breakpointSegments[i]->getBoundedLength() > minUpdateLength) {
// 			breakpointSegments[i]->updateBounds();
// 		}
		breakpointSegments[i]->updateBounds(minUpdateLength);
	}
}

void BreakpointGraph::resetScores() {
	std::for_each(breakpointSegments.begin(), breakpointSegments.end(),
				  std::mem_fun(&BreakpointSegment::resetScores));
}

void BreakpointGraph::setBestBreakpoints() {
	std::for_each(breakpointSegments.begin(), breakpointSegments.end(),
				  std::mem_fun(&BreakpointSegment::setBestBreakpoint));
}

void BreakpointGraph::refineMap() {
	std::for_each(breakpointSegments.begin(), breakpointSegments.end(),
				  std::mem_fun(&BreakpointSegment::split));
}

size_t BreakpointGraph::maxSegmentBoundedLength() {
	size_t maxLength = 0;
	for (std::vector<BreakpointSegment*>::const_iterator it = breakpointSegments.begin();
		 it != breakpointSegments.end(); ++it) {
		maxLength = std::max(maxLength, (*it)->getBoundedLength());
	}
	return maxLength;
}

void BreakpointGraph::updateSegmentScores() {
	VertexColorMap vertexColors = get(vertex_color, g);
	EdgeColorMap edgeColors = get(edge_color, g);
	Score totalScore = 0;
	undirected_dfs(g, ScoringVisitor(g), vertexColors, edgeColors);
	undirected_dfs(g, TotalScoreVisitor(g, totalScore), vertexColors, edgeColors);
	undirected_dfs(g, TracebackVisitor(g), vertexColors, edgeColors);

	std::cerr << "Total score: " << totalScore << '\n';
}

void BreakpointGraph::calcMidpointScore(util::Matrix<double>& distances) {
	BreakpointSegment::setResolution(3);
	initBounds();
	resetScores();
	readAlignments(distances);

	Score total = 0;
	
	AlignedSegmentsMap alignedSegmentsMap = get(aligned_segments_t(), g);
	EdgeIndexMap edgeIndexMap = get(edge_index_t(), g);	
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		AlignedSegments* a = alignedSegmentsMap[*ei];
		size_t midPos1 = a->getSegment1()->getNumIndices() / 2;
		size_t midPos2 = a->getSegment2()->getNumIndices() / 2;
		total += a->score(midPos1, midPos2);
	}

	std::cerr << "Midpoint break score: " << total << '\n';
}

void BreakpointGraph::calcBestBreakpoints(util::Matrix<double>& distances,
										  size_t maxSegmentPositions) {
	BreakpointSegment::setResolution(maxSegmentPositions);
	initBounds();
	resetScores();
	
	for (size_t iteration = 1; /* no check */; ++iteration) {
		std::cerr << "Iteration " << iteration << '\n';
		size_t maxLength = maxSegmentBoundedLength();
		std::cerr << "Max length = " << maxLength << '\n';

		size_t nextMaxLength = 2 * maxLength / (maxSegmentPositions - 1);
		
		readAlignments(distances);

		std::cerr << "Updating scores...\n";
		updateSegmentScores();

		if (maxLength + 1 <= maxSegmentPositions) {
			break;
		}

		std::cerr << "Updating bound...\n";
		updateBounds(std::max(nextMaxLength, maxSegmentPositions - 1));
		resetScores();
	}

	setBestBreakpoints();
}

void BreakpointGraph::writeMap(std::ostream& strm) {
	typedef std::vector<SegmentSet*>::const_iterator SegmentSetIt;
	for (SegmentSetIt it = segmentSets.begin(); it != segmentSets.end(); ++it) {
		strm << **it;
	}
}

void BreakpointGraph::writeBreakpointSegments(std::ostream& strm) {
	for (size_t i = 0; i < breakpointSegments.size(); ++i) {
		BreakpointSegment* seg = breakpointSegments[i];
		strm << seg->getNum() << '\t'
			 << seg->getGenome()->getName() << '\t'
			 << seg->getChrom() << '\t'
			 << seg->getStart() << '\t'
			 << seg->getEnd() << '\n';
	}
}
	
void BreakpointGraph::writeBreakpoints(std::ostream& strm) {
	for (size_t i = 0; i < breakpointSegments.size(); ++i) {
		BreakpointSegment* seg = breakpointSegments[i];
		strm << seg->getNum() << '\t'
			 << seg->getRelBreakpoint() << '\n';
	}
}

void BreakpointGraph::writeAlignedSegments(std::ostream& stream) {
	typedef std::vector<SegmentSet*>::const_iterator SegmentSetIt;
	for (SegmentSetIt it = segmentSets.begin(); it != segmentSets.end(); ++it) {
		(*it)->writeAlignedSegments(stream);
	}
}
