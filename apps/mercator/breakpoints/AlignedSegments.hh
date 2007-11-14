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

#ifndef __ALIGNED_SEGMENTS_HH__
#define __ALIGNED_SEGMENTS_HH__

#include "types.hh"

#include "bio/alignment/MultipleAlignment.hh"
using bio::alignment::MultipleAlignment;
#include "bio/alignment/AmbiguousDNAScoringMatrix.hh"
using bio::alignment::AmbiguousDNAScoringMatrix;
#include "util/matrix.hh"
using util::Matrix;
#include "filesystem.hh"

#include "BreakpointSegment.hh"

class AlignedSegments {
private:
	static filesystem::Path alignmentsDir;
	static std::string alignmentFilename;

	static const AmbiguousDNAScoringMatrix<Score>* matrix;
	static Score spaceScore;
	static Score gapScore;
	
	BreakpointSegment* segment1;
	BreakpointSegment* segment2;
	bool prefix1;
	bool prefix2;
	bool colinear;

	std::vector<size_t> pos1FromIndex2;
	std::vector<size_t> pos2FromIndex1;
	std::vector<Score> scoreFromIndex1;
	std::vector<Score> scoreFromIndex2;

	Score totalScore;
	
	bool bad_alignment;
	
	void processAlignment(MultipleAlignment& alignment);

	size_t getPosition1(size_t index1) const;
	size_t getPosition2(size_t index2) const;
	
public:
	static void setAlignmentsDir(const std::string& dir) {
		alignmentsDir = dir;
	}
	
	static void setAlignmentFilename(const std::string& filename) {
		alignmentFilename = filename;
	}
	
	static void setScores(const AmbiguousDNAScoringMatrix<Score>& matrix,
						  Score spaceScore, Score gapScore);
	
	AlignedSegments(BreakpointSegment* segment1,
					BreakpointSegment* segment2,
					bool prefix1,
					bool prefix2);

	BreakpointSegment* getSegment1() const { return segment1; }
	BreakpointSegment* getSegment2() const { return segment2; }
	bool isPrefix1() const { return prefix1; }
	bool isPrefix2() const { return prefix2; }
	bool isColinear() const { return colinear; }
	
	// Returns true if THIS and OTHER align the same two segments,
	// possibly in opposite orientations
	bool hasSameSegments(const AlignedSegments& other) const;
	
	// Reads the pairwise alignment of the two segments
	void readAlignment(size_t alignmentNum, Matrix<double>& distances);

	// Returns the score of this alignment if segment1 is cut at POS1
	// and segment2 is cut at POS2
	Score score(size_t pos1, size_t pos2) const;
};

inline size_t
AlignedSegments::getPosition1(size_t index1) const {
	return (prefix1 ?
			segment1->getPosition(index1) :
			segment1->getLength() - segment1->getPosition(index1));
}

inline size_t
AlignedSegments::getPosition2(size_t index2) const {
	return (prefix2 ?
			segment2->getPosition(index2) :
			segment2->getLength() - segment2->getPosition(index2));
}

#endif // __ALIGNED_SEGMENTS_HH__
