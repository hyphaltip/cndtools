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
#include "bio/agp/AGPBackwardMapper.hh"

namespace bio { namespace agp {

	struct AGPRecComparer {
		bool operator()(const formats::agp::Record& x,
						const formats::agp::Record& y) const {
			int chromComp = x.getChrom().compare(y.getChrom());
			return chromComp < 0
				or (chromComp == 0
					and (x.getContigStart() < y.getContigStart()
						 or (x.getContigStart()  == y.getContigStart()
							 and x.getContigEnd() < y.getContigEnd())));
		}
	};

	struct AGPIntervalComparer {
		bool operator()(const formats::agp::Record& x,
						const genome::Interval& y) const{
			int chromComp = x.getChrom().compare(y.getChrom());
			return chromComp < 0
				or (chromComp == 0 && x.getContigEnd() <= y.getStart());
		}
	};
	
	AGPBackwardMapper::AGPBackwardMapper(std::istream& strm) {
		// Read in AGP file
		formats::agp::InputStream agpStream(strm);
		formats::agp::Record agpRec;
		while (agpStream >> agpRec) {
			if (!agpRec.isGap()) {
				agpRecs.push_back(agpRec);
			}
		}

		// Sort AGP Records by chrom interval
		std::sort(agpRecs.begin(), agpRecs.end(), AGPRecComparer());
	}
	
	void AGPBackwardMapper::map(const genome::Interval& i,
								std::vector<genome::BasicInterval>& mapped) {
		genome::BasicInterval query(i);
		std::vector<formats::agp::Record>::const_iterator start
			= std::lower_bound(agpRecs.begin(), agpRecs.end(),
							   query, AGPIntervalComparer());
		if (start == agpRecs.end()) {
			return;
		}

		// Find last AGP record that this interval overlaps
		std::vector<formats::agp::Record>::const_iterator end = start;
		while (end != agpRecs.end() and overlap(*end, i)) {
			++end;
		}			 

		std::vector<formats::agp::Record>::const_iterator it;
		for (it = start; it != end; ++it) {
			genome::BasicInterval untransformedInt(i);

			genome::BasicInterval contig_interval = it->getTargetInterval();
			genome::BasicInterval source_interval = it->getSourceInterval();
			
			int startOffset = i.getStart() < contig_interval.getStart() ? 
				0 :
				i.getStart() - contig_interval.getStart();
			int endOffset = i.getEnd() < contig_interval.getEnd() ?
				i.getEnd() - contig_interval.getStart() :
				contig_interval.getLength();

			untransformedInt.setChrom(it->getSourceAccession());
			
			if (it->getSourceOrientation() == '+') {
				untransformedInt.setStart(source_interval.getStart() + startOffset);
				untransformedInt.setEnd(source_interval.getStart() + endOffset);
			} else {
				untransformedInt.setStart(source_interval.getEnd() - endOffset);
				untransformedInt.setEnd(source_interval.getEnd() - startOffset);
				untransformedInt.setStrand(untransformedInt.getStrand().opposite());
			}

			mapped.push_back(untransformedInt);
		}
		
	}

	bool AGPBackwardMapper::overlap(const formats::agp::Record& rec,
									const genome::Interval& i) {
		return rec.getContigEnd() > i.getStart()
			and rec.getContigStart() <= i.getEnd()
			and rec.getChrom() == i.getChrom();
	}
		
} }
