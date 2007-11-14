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

#ifndef __CHROMOSOME_HH__
#define __CHROMOSOME_HH__

#include "types.hh"

class Chromosome {
public:
	Chromosome(Genome* genome,
			   const string& name,
			   const GenomicDist length);

	virtual ~Chromosome() {}

	Genome* getGenome() const;
	const string& getName() const;
	string getUniqueName() const;
	GenomicDist getLength() const;
	GenomicDist getGenomeStart() const;
	GenomicDist getCoverage() const;
	float getPctCoverage() const;
	bool isReversed() const;
	bool isPartOfAssembled() const;
	
	void setGenomeStart(GenomicDist start);

	void addAnchor(Anchor* a);
	void initAnchors();
	
	size_t getNumAnchors() const;
	size_t getNumAnchorsInCliques() const;
	size_t getNumAnchorsInRuns() const;
	size_t getNumCliques() const;
	size_t getNumRuns() const;

	Anchor* getFirstAnchor() const;
	Anchor* getLastAnchor() const;
	Clique* getFirstClique() const;
	Clique* getLastClique() const;
	Run* getFirstRun() const;
	Run* getLastRun() const;

	// Reverse this chromosome, changing the coordinates and signs of
	// the anchors on it
	void reverse();

	// Unflip all anchors on this chromosome
	void unflipAnchors();
	
	// Returns the number of anchors marked as repetitive in this
	// chromosome
	size_t getNumAnchorsRepetitive() const;

	// Write the permutation of run numbers on this chromosome to STRM
	void writeRunPerm(std::ostream& stream) const;

	// Write AGP lines for this chromosome
	virtual void writeAGPLines(std::ostream& stream, size_t& recNum) const;
	
protected:
	vector<Anchor*> anchors;
	Genome* genome;
	string name;
	size_t num;
	GenomicDist length;
	GenomicDist genomeStart;
	bool reversed;
};

struct ChromLengthSorter {
	bool operator()(const Chromosome* c1, const Chromosome* c2) const {
		return c1->getLength() < c2->getLength();
	}
};

#endif // __CHROMOSOME_HH__
