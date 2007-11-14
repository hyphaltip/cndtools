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

#include "assembled.hh"
#include "anchor.hh"

#include "bio/formats/agp/Record.hh"

GenomicDist Assembled::padding = 100;

void Assembled::setPadding(GenomicDist n) {
	padding = n;
}

GenomicDist Assembled::getPadding() {
	return padding;
}

Assembled::Assembled(Genome* genome,
					 const string& name) :
	Chromosome(genome, name, 0)
{}

void Assembled::addChrom(const Chromosome* chrom) {
	// Record name and orientation of chromosome
	chromNames.push_back(chrom->getName());
	chromOrients.push_back((chrom->isReversed() ? '-' : '+'));
	chromLengths.push_back(chrom->getLength());

	// Copy the anchors into this chromosome

	// Calculate the offset for the anchors to be added depending on whether
	// this is the first chromosome added or not (affects padding)
	GenomicDist offset = (chromNames.size() == 1 ? 0 : length + padding);
	
	for (Anchor* a = chrom->getFirstAnchor(); a != NULL; a = a->nextAnchor()) {
		a->setChrom(this);
		a->shift(offset);
		if (!anchors.empty()) {
			a->setPrev(anchors.back());
			anchors.back()->setNext(a);
		}
		anchors.push_back(a);
	}

	// Update the length of the assembled chromosome
	length += (chromNames.size() == 1 ? 0 : padding) + chrom->getLength();
}

size_t Assembled::numChroms() const {
	return chromNames.size();
}

void Assembled::writeAGPLines(ostream& stream, size_t& recNum) const {
	GenomicDist startCoord = 1;
	for (size_t i = 0; i < numChroms(); ++i) {
		
		if (i != 0 && padding != 0) {
			stream << bio::formats::agp::Record(name,
												"",
												startCoord,
												startCoord + padding - 1,
												recNum,
												'N',
												padding,
												"contig",
												true);
			++recNum;
			startCoord += padding;
		}
		
		stream << bio::formats::agp::Record(name,
											"",
											startCoord,
											startCoord + chromLengths[i] - 1,
											recNum,
											'D',
											chromNames[i],
											1,
											chromLengths[i],
											chromOrients[i]);
		++recNum;
		startCoord += chromLengths[i];
	}
}
