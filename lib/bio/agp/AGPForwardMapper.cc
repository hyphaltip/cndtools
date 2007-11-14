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

#include "bio/formats/agp/InputStream.hh"
#include "bio/agp/AGPForwardMapper.hh"

namespace bio { namespace agp {
	
	AGPForwardMapper::AGPForwardMapper(std::istream& stream) {
		// Read in AGP file
		formats::agp::InputStream agpStream(stream);
		formats::agp::Record agpRec;
		while (agpStream >> agpRec) {
			if (!agpRec.isGap()) {
				contigMap[agpRec.getSourceAccession()].push_back(agpRec);
			}
			chromSizes[agpRec.getChrom()] += agpRec.getContigLength();
		}
	}
	
	void AGPForwardMapper::map(const genome::Interval& i,
							   std::vector<genome::BasicInterval>& mapped) {
		AGPMap::iterator it = contigMap.find(i.getChrom());
		if (it == contigMap.end()) {
			return;
		}
		
		AGPList& mappingList = it->second;

		for (AGPList::iterator mapping = mappingList.begin();
			 mapping != mappingList.end(); ++mapping) {
			if (not overlap(i, *mapping)) {
				continue;
			}
			
			genome::BasicInterval j;
			
			j.setChrom(mapping->getChrom());
			
			if (mapping->getSourceOrientation() == '+') {
				if (i.getStart() > mapping->getSourceStart() - 1) {
					j.setStart(i.getStart()
							   + mapping->getContigStart()
							   - mapping->getSourceStart());
				} else {
					j.setStart(mapping->getContigStart() - 1);
				}
				if (i.getEnd() < mapping->getSourceEnd()) {
					j.setEnd(i.getEnd()
							 + mapping->getContigStart()
							 - mapping->getSourceStart());
				} else {
					j.setEnd(mapping->getContigEnd());
				}
				j.setStrand(i.getStrand());
			} else {
				if (i.getEnd() < mapping->getSourceEnd()) {
					j.setStart(mapping->getSourceEnd()
							   - i.getEnd()
							   + mapping->getContigStart() - 1);
				} else {
					j.setStart(mapping->getContigStart() - 1);
				}
				if (i.getStart() > mapping->getSourceStart() - 1) {
					j.setEnd(mapping->getSourceEnd()
							 - i.getStart()
							 + mapping->getContigStart() - 1);
				} else {
					j.setEnd(mapping->getContigEnd());
				}
				j.setStrand(i.getStrand().opposite());
			}
			
			mapped.push_back(j);
		}
	}

	genome::Distance
	AGPForwardMapper::getChromSize(const std::string& chrom) const {
		LengthMap::const_iterator it = chromSizes.find(chrom);
		if (it == chromSizes.end()) {
			throw std::runtime_error("Chromosome not in AGP file: " + chrom);
		}
		return it->second;
	}
	
	bool AGPForwardMapper::overlap(const genome::Interval& i,
								   const formats::agp::Record& rec) {
		return i.getStart() < rec.getSourceEnd()
			and i.getEnd() >= rec.getSourceStart()
			and i.getChrom() == rec.getSourceAccession();
	}
	
} }
