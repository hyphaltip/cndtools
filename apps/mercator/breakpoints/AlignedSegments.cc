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

#include "AlignedSegments.hh"
#include "BreakpointSegment.hh"
#include "Genome.hh"
#include "Segment.hh"

#include "bio/alphabet/Nucleotide.hh"
#include "bio/alignment/BasicNamedMultipleAlignment.hh"
#include "bio/formats/fasta/InputStream.hh"
#include "util/string.hh"
#include "util/stl.hh"
using bio::alphabet::Nucleotide;
using bio::alignment::BasicNamedMultipleAlignment;
using util::stl::print_elements;
using namespace bio;

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>

filesystem::Path AlignedSegments::alignmentsDir;
std::string AlignedSegments::alignmentFilename;

const AmbiguousDNAScoringMatrix<Score>* AlignedSegments::matrix = NULL;
Score AlignedSegments::spaceScore = 0;
Score AlignedSegments::gapScore = 0;

void AlignedSegments::setScores(const AmbiguousDNAScoringMatrix<Score>& matrix,
								Score spaceScore, Score gapScore) {
	AlignedSegments::matrix = &matrix;
	AlignedSegments::spaceScore = spaceScore;
	AlignedSegments::gapScore = gapScore;
}

AlignedSegments::AlignedSegments(BreakpointSegment* segment1,
								 BreakpointSegment* segment2,
								 bool prefix1,
								 bool prefix2)
	: segment1(segment1),
	  segment2(segment2),
	  prefix1(prefix1),
	  prefix2(prefix2),
	  colinear(segment1->isColinearWith(*segment2))
{
}

void AlignedSegments::readAlignment(size_t alignmentNum,
									Matrix<double>& distances) {
	// Check if either breakpoint segment has changed in size.  If
	// not, we do not need to realign these segments
	if (not (segment1->isChanged() or segment2->isChanged())) {
		return;
	}
	
 	std::cerr << "Reading alignment:"
 			  << '\t' << alignmentNum
 			  << '\t' << *segment1
 			  << '\t' << *segment2
 			  << '\n';

	BasicNamedMultipleAlignment alignment;
	std::string num = util::string::toString(alignmentNum);
	filesystem::Path filename(alignmentsDir / num / alignmentFilename);
	filesystem::InputFileStream alignmentFile;	

	try {
		alignmentFile.open(filename);
	} catch (const std::runtime_error& e) {
		// Do nothing
	}

	formats::fasta::InputStream fastaStream(alignmentFile);
	fastaStream >> alignment;

	if (alignment.getNumSeqs() != 2 or
		(alignment.getNumSeqs() == 2 and
		 (alignment.getSeqLen(0) != segment1->getLength() or
		  alignment.getSeqLen(1) != segment2->getLength()))) {
		std::cerr << "Warning: Invalid alignment file: "
				  << filename.toString()
				  << '\n';
// 		alignment.clear();
// 		std::string seq1 = segment1->getSeq((prefix1 ? '+' : '-')) +
// 			std::string(segment2->getLength(), '-');
// 		std::string seq2 = std::string(segment1->getLength(), '-') +
// 			segment2->getSeq((prefix2 ? '+' : '-'));
// 		alignment.addSeq(seq1, "seq1");
// 		alignment.addSeq(seq2, "seq2");
		bad_alignment = true;
	} else {
		bad_alignment = false;
		processAlignment(alignment);
	}
}

void AlignedSegments::processAlignment(MultipleAlignment& alignment) {	
	std::vector<size_t> pos1FromPos2(segment2->getLength() + 1);
	std::vector<size_t> pos2FromPos1(segment1->getLength() + 1);
	std::vector<Score> scoreFromPos1(segment1->getLength() + 1);
	std::vector<Score> scoreFromPos2(segment2->getLength() + 1);
	
	pos2FromIndex1.resize(segment1->getNumIndices());
	scoreFromIndex1.resize(segment1->getNumIndices());
	pos1FromIndex2.resize(segment2->getNumIndices());
	scoreFromIndex2.resize(segment2->getNumIndices());

	pos1FromPos2.at(0) = 0;
	pos2FromPos1.at(0) = 0;
	scoreFromPos1.at(0) = 0;
	scoreFromPos2.at(0) = 0;
	
	bool inGap1 = false;
	bool inGap2 = false;
	size_t firstPrefix = 0;
	size_t secondPrefix = 0;
	Score partialScore = 0;

// 	std::cerr << alignment.getSeq(0) << '\n';
// 	std::cerr << alignment.getSeq(1) << '\n';
	
	for (size_t i = 0; i < alignment.getNumCols(); ++i) {
		char c1 = alignment.getChar(0, i);
		char c2 = alignment.getChar(1, i);

		firstPrefix += (c1 == '-' ? 0 : 1);
		secondPrefix += (c2 == '-' ? 0 : 1);

		if (c1 != '-' and c2 != '-') {
			partialScore += matrix->getCharScore(c1, c2);
		} else if (c1 == '-') {
			partialScore += (inGap1 ? spaceScore : gapScore + spaceScore);
		} else if (c2 == '-') {
			partialScore += (inGap2 ? spaceScore : gapScore + spaceScore);
		} else {
			throw std::runtime_error("Alignment has gap aligned to gap");
		}

		inGap1 = (c1 == '-');
		inGap2 = (c2 == '-');

		if (not inGap1) {
			scoreFromPos1.at(firstPrefix) = partialScore;
			pos2FromPos1.at(firstPrefix) = secondPrefix;
		}
		if (not inGap2) {
			scoreFromPos2.at(secondPrefix) = partialScore;
			pos1FromPos2.at(secondPrefix) = firstPrefix;
		}

// 		std::cerr << c1 << '\t' << c2 << '\t'
// 				  << firstPrefix << '\t' << secondPrefix << '\t'
// 				  << partialScore << '\n';
	}
	totalScore = partialScore;

	for (size_t index1 = 0; index1 < segment1->getNumIndices(); ++index1) {
		size_t pos = getPosition1(index1);
		pos2FromIndex1.at(index1) = pos2FromPos1.at(pos);
		scoreFromIndex1.at(index1) = scoreFromPos1.at(pos);
	}
	for (size_t index2 = 0; index2 < segment2->getNumIndices(); ++index2) {
		size_t pos = getPosition2(index2);
		pos1FromIndex2.at(index2) = pos1FromPos2.at(pos);
		scoreFromIndex2.at(index2) = scoreFromPos2.at(pos);
	}

//  	std::cerr << "AlignedSegment:" << *segment1 << " " << *segment2 << '\n';
	
// 	std::cerr << "pos1FromPos2\n";
// 	print_elements(std::cerr, pos1FromPos2);
// 	std::cerr << "pos2FromPos1\n";
// 	print_elements(std::cerr, pos2FromPos1);
// 	std::cerr << "scoreFromPos1\n";
// 	print_elements(std::cerr, scoreFromPos1);
// 	std::cerr << "scoreFromPos2\n";
// 	print_elements(std::cerr, scoreFromPos2);

// 	std::cerr << "pos1FromIndex2\n";
// 	print_elements(std::cerr, pos1FromIndex2);
// 	std::cerr << "pos2FromIndex1\n";
// 	print_elements(std::cerr, pos2FromIndex1);
// 	std::cerr << "scoreFromIndex1\n";
// 	print_elements(std::cerr, scoreFromIndex1);
// 	std::cerr << "scoreFromIndex2\n";
// 	print_elements(std::cerr, scoreFromIndex2);
}

Score AlignedSegments::score(size_t index1, size_t index2) const {
	try {
		if (bad_alignment) { return 0; }
		
		size_t pos1 = getPosition1(index1);
		size_t pos2 = getPosition2(index2);
		
		size_t prefix1 = pos1FromIndex2.at(index2);
		size_t prefix2 = pos2FromIndex1.at(index1);

		if (colinear) {
			if (pos2 > prefix2) {
				return scoreFromIndex1.at(index1)
					+ spaceScore * (pos2 - prefix2)
					+ (totalScore - scoreFromIndex2.at(index2))
					+ spaceScore * (prefix1 - pos1);
			} else if (pos1 >= prefix1) {
				return scoreFromIndex2.at(index2) 
					+ spaceScore * (pos1 - prefix1)
					+ (totalScore - scoreFromIndex1.at(index1))
					+ spaceScore * (prefix2 - pos2);
			} else {
				throw std::runtime_error("pos1 < prefix1");
			}
		} else {
			if (pos2 > prefix2) {
				return scoreFromIndex1.at(index1) + spaceScore * (pos2 - prefix2);
			} else if (pos1 >= prefix1) {
				return scoreFromIndex2.at(index2) + spaceScore * (pos1 - prefix1);
			} else {
				throw std::runtime_error("pos1 < prefix1");
			}
		}
	} catch (std::runtime_error& e) {
		std::cerr << "Error: score(): " << e.what() << '\n';
		return 0;
	}
}
