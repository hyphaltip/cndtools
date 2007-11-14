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
#include <limits>

#include "util/string.hh"
#include "filesystem.hh"
#include "bio/homologymap/HomologyMapper.hh"
#include "bio/formats/fasta/InputStream.hh"
using util::string::toString;
using namespace filesystem;

namespace bio { namespace homologymap {
	
	const std::string HomologyMapper::GENOMES_FILENAME = "genomes";
	const std::string HomologyMapper::MAP_FILENAME = "map";
	
	HomologyMapper::HomologyMapper(const filesystem::Path& alignDir,
								   const std::string& sourceGenome,
								   const std::string& targetGenome)
		: alignDir(alignDir),
		  lastSegmentNum(std::numeric_limits<size_t>::max()),
		  sourceGenome(sourceGenome),
		  targetGenome(targetGenome) {

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
		size_t index = 0;
		while (genomeFile >> genome) {
			genomeIndices[genome] = index;
			++index;
		}
		
		// Read map
		std::cerr << "Reading map...\n";
		hmap.read(mapFile);
	}
	
	void HomologyMapper::map(const genome::Interval& i,
							 std::vector<genome::BasicInterval>& mapped) {
		map(i, mapped, sourceGenome, targetGenome);
	}
	
	size_t HomologyMapper::getIndex(const std::string& g) const {
		util::stl::hash_map<std::string, size_t>::const_iterator it;
		it = genomeIndices.find(g);
		if (it == genomeIndices.end()) {
			throw std::runtime_error(g + " not in map");
		}
		return it->second;
	}
			
	bool HomologyMapper::readAlignment(size_t segmentNum) {
		if (segmentNum == lastSegmentNum) { return true; }
		
		// Read in alignment file
		filesystem::InputFileStream segFile;
		try {
			segFile.open(alignDir / toString(segmentNum) / "mavid.mfa");
		} catch (const std::runtime_error& e) {
			return false;
		}
		formats::fasta::InputStream fastaStream(segFile);
		fastaStream >> align;
		lastSegmentNum = segmentNum;
		return true;
	}

	void HomologyMapper::flipInterval(alignment::Interval& i, size_t length) {
		std::swap(i.start, i.end);
		i.start = length - i.start;
		i.end = length - i.end;
	}

	void HomologyMapper::getInterval(genome::Interval* segSourceInt,
									 genome::Interval* segTargetInt,
									 int alignSourceNum,
									 int alignTargetNum,
									 const genome::Interval& sourceInterval,
									 genome::BasicInterval& targetInterval) {
		alignment::Interval alignSourceInt;
		if (segSourceInt->contains(sourceInterval.getStartCoord())) {
			alignSourceInt.start = sourceInterval.getStart() - segSourceInt->getStart();
		} else {
			alignSourceInt.start = 0;
		}
		
		if (segSourceInt->contains(sourceInterval.getEndCoord())) {
			alignSourceInt.end = sourceInterval.getEnd() - segSourceInt->getStart();
		} else {
			alignSourceInt.end = segSourceInt->getLength();
		}
				
		// Flip source interval if it is on the reverse strand
		if (segSourceInt->getStrand().isReverse()) {
			flipInterval(alignSourceInt, segSourceInt->getLength());
		}
		
		// Get column interval in alignment corresponding to
		// this source interval
		alignment::Interval columnInt =
			align.getColumnInterval(alignSourceNum,
									alignSourceInt.start,
									alignSourceInt.end);
				
		// Get sequence interval in target
		alignment::Interval alignTargetInt =
			align.getSeqInterval(alignTargetNum,
								 columnInt.start,
								 columnInt.end);
						
		// Flip target interval if it is on the reverse strand
		if (segTargetInt->getStrand().isReverse()) {
			flipInterval(alignTargetInt, segTargetInt->getLength());
		}

		targetInterval.setChrom(segTargetInt->getChrom());
		targetInterval.setStart(segTargetInt->getStart() + alignTargetInt.start);
		targetInterval.setEnd(segTargetInt->getStart() + alignTargetInt.end);
		targetInterval.setStrand(segSourceInt->getStrand()
								 == segTargetInt->getStrand()
								 ? sourceInterval.getStrand()
								 : sourceInterval.getStrand().opposite());
	}
	
	void HomologyMapper::map(const genome::Interval& sourceInterval,
							 std::vector<genome::BasicInterval>& mapped,
							 const std::string& source,
							 const std::string& target) {
		size_t sourceIndex = getIndex(source);
		size_t targetIndex = getIndex(target);
		
		std::vector<homologymap::Segment*> segments;
		hmap.getSegments(sourceIndex, sourceInterval, segments);

		if (segments.empty()) { return; }

		genome::BasicInterval targetInterval;
		for (size_t i = 0; i < segments.size(); ++i) {
			Segment* seg = segments[i];
			if (not seg->hasGenome(targetIndex)) { continue; }

			genome::Interval* segSourceInt = seg->intervals[sourceIndex];
			genome::Interval* segTargetInt = seg->intervals[targetIndex];
			
			if (not readAlignment(seg->num)) {
				std::cerr << "Warning: Could not open alignment file for "
						  << "segment " << seg->num << '\n';
				continue;
			}
			
			int alignSourceNum = align.getSeqNum(source);
			int alignTargetNum = align.getSeqNum(target);

			getInterval(segSourceInt, segTargetInt,
						alignSourceNum, alignTargetNum,
						sourceInterval, targetInterval);
			
			if (targetInterval.getLength() > 0) {
				mapped.push_back(targetInterval);
			}
		}
	}		
		
	void HomologyMapper::setTargetGenome(const std::string& targetGenome) {
		this->targetGenome = targetGenome;
	}

	void HomologyMapper::setSourceGenome(const std::string& sourceGenome) {
		this->sourceGenome = sourceGenome;
	}

} }
