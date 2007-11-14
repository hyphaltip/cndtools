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

#ifndef __ASSEMBLED_HH__
#define __ASSEMBLED_HH__

#include "types.hh"
#include "chromosome.hh"

class Assembled : public Chromosome {
public:
	// Set the padding between original chromosomes in this assembly
	static void setPadding(GenomicDist n);
	
	// Get the padding between original chromosomes in this assembly
	static GenomicDist getPadding();

	// Make a new assembled chromosome in genome GENOME with name NAME
	Assembled(Genome* genome, const string& name);

	// Add chromosome CHROM as part of this assembled chromosome.
	// Modifies the coordinates of the anchors in CHROM to reflect
	// their new assembled coordinates and sets the chromosome of the
	// anchors to be THIS
	void addChrom(const Chromosome* chrom);

	// Returns the number of original chromosomes that make up this
	// assembled chromosome
	size_t numChroms() const;

	// Output the AGP records that can be used to assembled this
	// chromosome onto the stream STRM, with record numbers starting
	// at RECNUM.  RECNUM will incremented by the number of records
	// written upon returning.
	void writeAGPLines(ostream& stream, size_t& recNum) const;

private:
	// The number of N's to put in between original chromosome sequences 
	static GenomicDist padding;
	
	// The names of the original chromosomes that make up this
	// assembled chromosome
	vector<string> chromNames;

	// The orientations of the original chromosomes that make up this
	// assembled chromosome (+ or -)
	vector<char> chromOrients;
	
	// The lengths of the original chromosomes that make up this
	// assembled chromosome
	vector<GenomicDist> chromLengths;
};

#endif // __ASSEMBLED_HH__
