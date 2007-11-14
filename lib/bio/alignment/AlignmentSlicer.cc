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

#include <iostream>

#include "bio/alignment/AlignmentSlicer.hh"
#include "bio/alphabet/AmbiguousNucleotide.hh"
#include "bio/genome/BasicInterval.hh"
#include "bio/formats/fasta.hh"
#include "util/string.hh"
using util::string::toString;
using namespace filesystem;
using namespace bio::formats;

namespace bio { namespace alignment {

	const std::string AlignmentSlicer::GENOMES_FILENAME = "genomes";
	const std::string AlignmentSlicer::MAP_FILENAME = "map";
	
	AlignmentSlicer::
	AlignmentSlicer(const Path& alignDir,
					const std::string& sourceGenome,
					char noMapChar)
		: alignDir(alignDir),
		  sourceGenome(sourceGenome),
		  noMapChar(noMapChar),
		  genomes(),
		  indexMap(),
		  sourceIndex(0),
		  map(),
		  align(),
		  lastSegNum(std::numeric_limits<size_t>::max())
	{
		// Open genomes and map files
		InputFileStream genomeFile(alignDir / GENOMES_FILENAME);
		InputFileStream mapFile(alignDir / MAP_FILENAME);

		// Check for existence of alignment dir
		if (not alignDir.exists()) {
			throw std::runtime_error("Alignment directory does not exist: " +
									 alignDir.toString());
		}

		// Read in genomes
		std::string genome;
		while (genomeFile >> genome) {
			indexMap[genome] = genomes.size();
			genomes.push_back(genome);
		}

		// Set source genome index
		if (indexMap.find(sourceGenome) == indexMap.end()) {
			throw std::runtime_error(sourceGenome + " not in map");
		}
		sourceIndex = indexMap[sourceGenome];

		// Read map
		std::cerr << "Reading map...\n";
		map.read(mapFile);
	}

	bool
	AlignmentSlicer::
	readAlignment(size_t segNum) {
		if (segNum == lastSegNum) {
			return align.getNumSeqs() > 0;
		}

		std::cerr << "Reading alignment...\n";
		lastSegNum = segNum;
		align.clear();
		InputFileStream segFile;
		try {
			segFile.open(alignDir / toString(segNum) / "mavid.mfa");
		} catch (const std::runtime_error& e) {
			std::cerr << "Warning: Could not open alignment file for "
					  << "segment " << segNum << '\n';
			return false;
		}
		fasta::InputStream fastaStream(segFile);
		fastaStream >> align;
		return true;
	}

	BasicNamedMultipleAlignment
	AlignmentSlicer::
	getSlice(const bio::genome::Interval& sourceInterval) {
		typedef std::vector<bio::genome::BasicInterval> IntervalList;
		typedef std::vector<std::string> StringList;
		std::vector<IntervalList> intervals(genomes.size());
		std::vector<StringList> slices(genomes.size());
		
		std::vector<bio::homologymap::Segment*> segments;
		
		//std::cerr << "sourceInterval: " << sourceInterval << '\n';

		segments.clear();
		map.getSegments(sourceIndex, sourceInterval, segments);
		if (segments.empty()) {
			throw std::runtime_error("Interval not in map: " +
									 toString(sourceInterval));
		}

		if ((sourceInterval.getStartCoord() <
			 segments.front()->intervals[sourceIndex]->getStartCoord()) or
			(segments.back()->intervals[sourceIndex]->getEndCoord() <
			 sourceInterval.getEndCoord())) {
			std::cerr << "Warning: Interval is only partially in map: "
					  << sourceInterval << '\n';
		}

		for (size_t i = 0; i < segments.size(); ++i) {
			if (not readAlignment(segments[i]->num)) {
				continue;
			}
			
			int alignSourceNum = align.getSeqNum(sourceGenome);

			bio::genome::Interval* segSourceInt =
				segments[i]->intervals[sourceIndex];

			//std::cerr << "segSourceInt: " << *segSourceInt << '\n';
				
			bio::alignment::Interval alignSourceInt;
			if (i == 0 and
				segSourceInt->contains(sourceInterval.getStartCoord())) {
				alignSourceInt.start = sourceInterval.getStart() - segSourceInt->getStart();
			} else {
				alignSourceInt.start = 0;
			}

			if (i == (segments.size() - 1) and
				segSourceInt->contains(sourceInterval.getEndCoord())) {
				alignSourceInt.end = sourceInterval.getEnd() - segSourceInt->getStart();
			} else {
				alignSourceInt.end = segSourceInt->getLength();
			}

// 			std::cerr << "sourceInt: " << alignSourceInt.start << "-"
// 			<< alignSourceInt.end << '\n';

			// Flip source interval if it is on the reverse strand
			if (segSourceInt->getStrand() == '-') {
				std::swap(alignSourceInt.start, alignSourceInt.end);
				alignSourceInt.start = segSourceInt->getLength() - alignSourceInt.start;
				alignSourceInt.end = segSourceInt->getLength() - alignSourceInt.end;
			}

// 			std::cerr << "sourceInt: " << alignSourceInt.start << "-"
// 					  << alignSourceInt.end << '\n';
				
			// Get column interval in alignment corresponding to
			// this source interval
			bio::alignment::Interval columnInt =
				align.getColumnInterval(alignSourceNum,
										alignSourceInt.start,
										alignSourceInt.end);

// 			std::cerr << "columnInt: " << columnInt.start << '-'
// 			<< columnInt.end << '\n';

			size_t numColumns = columnInt.end - columnInt.start;

			for (size_t targetIndex = 0; targetIndex < genomes.size();
				 ++targetIndex) {
				std::string targetGenome = genomes[targetIndex];
				int alignTargetNum = align.getSeqNum(targetGenome);
				bio::genome::Interval* segTargetInt = segments[i]->intervals[targetIndex];
				if (segTargetInt == NULL) {
					slices[targetIndex].push_back(std::string(numColumns, noMapChar));
					continue;
				}
				//std::cerr << "segTargetInt: " << *segTargetInt << '\n';
				
				// Get sequence interval in target
				bio::alignment::Interval alignTargetInt =
					align.getSeqInterval(alignTargetNum,
										 columnInt.start,
										 columnInt.end);

				//std::cerr << "targetInt: " << alignTargetInt.start << '-'
				//		  << alignTargetInt.end << '\n';

				// Skip if there is no sequence in target matching source
				if (alignTargetInt.start == alignTargetInt.end) {
					slices[targetIndex].push_back(std::string(numColumns, '-'));
					continue;
				}
				
				// Flip target interval if it is on the reverse strand
				if (segTargetInt->getStrand() == '-') {
					std::swap(alignTargetInt.start, alignTargetInt.end);
					alignTargetInt.start = segTargetInt->getLength() - alignTargetInt.start;
					alignTargetInt.end = segTargetInt->getLength() - alignTargetInt.end;
				}
				
				bio::genome::BasicInterval
					targetInt(segTargetInt->getChrom(),
							  segTargetInt->getStart() + alignTargetInt.start,
							  segTargetInt->getStart() + alignTargetInt.end,
							  segTargetInt->getStrand());

				//std::cerr << "targetInterval: " << targetInt << '\n';

				std::string targetSlice = align.getSubstring(alignTargetNum,
															 columnInt.start,
															 columnInt.end);

				if (sourceInterval.getStrand() != segSourceInt->getStrand()) {
					bio::alphabet::AmbiguousDNA.reverseComplementInPlace(targetSlice);
					targetInt.flip();
				}
				
				intervals[targetIndex].push_back(targetInt);
				slices[targetIndex].push_back(targetSlice);
			}
		}

		BasicNamedMultipleAlignment alignmentSlice;
		for (size_t i = 0; i < genomes.size(); ++i) {
			if (sourceInterval.getStrand().isReverse()) {
				std::reverse(intervals[i].begin(), intervals[i].end());
				std::reverse(slices[i].begin(), slices[i].end());
			}

			std::string name = genomes[i];
			for (size_t j = 0; j < intervals[i].size(); ++j) {
				bio::genome::Interval& interval = intervals[i][j];
				name += " "
					+ interval.getChrom() + ":"
					+ toString(interval.getStart()) + "-"
					+ toString(interval.getEnd())
					+ toString(interval.getStrand());
			}

			std::string sequence;
			for (size_t j = 0; j < slices[i].size(); ++j) {
				sequence += slices[i][j];
			}
			
			alignmentSlice.addSeq(sequence, name);
		}

		return alignmentSlice;
	}

} }
