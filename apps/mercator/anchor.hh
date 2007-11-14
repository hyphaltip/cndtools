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

#ifndef __ANCHOR_HH__
#define __ANCHOR_HH__

#include "types.hh"
#include "chromosome.hh"
#include "genome.hh"
#include "clique.hh"
#include "run.hh"

class Anchor {
public:

	Anchor(const string& name,
		   Chromosome* chrom,
		   const char strand,
		   const GenomicDist start,
		   const GenomicDist end,
		   const size_t isCoding);

	const string& getName() const;

	size_t getNameNum() const;

	void print() const;
	
	Chromosome* getChrom() const;
	Genome* getGenome() const;
	string getUniqueName() const;
	char getStrand() const;

	// Returns true if this anchor is on the + strand
	bool isForward() const;

	bool isMarked() const;
	void setMarked(const bool state);

	// Returns true if this anchor has (to any genome) two or more
	// hits with the maximum score, or at least REPEATNUM hits with
	// scores within REPEATPCT of the maximum hit score (considered
	// separately in each genome)
	bool isRepetitive(size_t repeatNum, float repeatPct) const;
	
	// Returns the distance from this anchor to the end of its
	// chromosome in the direction indicated by FORWARD
	GenomicDist getDistanceToChromEnd(const bool forward = true) const;
	
	size_t getGenomeNum() const;
	
	GenomicDist getStart() const;
	GenomicDist getEnd() const;
	GenomicDist getLength() const;
	GenomicDist getGenomeStart() const;
	GenomicDist getGenomeEnd() const;
	GenomicDist getGenomeMiddle() const;
	GenomicDist getDistanceTo(Anchor* other, bool forward) const;

	void setChrom(Chromosome* chrom);
	void setClique(Clique* c);
	void setNext(Anchor* a);
	void setPrev(Anchor* a);

	// Returns the number of edges incident to this anchor from GENOME
	int getNumEdges(const size_t genome) const;

	// Returns the maximum (over genomes) number of edges incident to
	// this anchor
	int getMaxEdges() const;
	
	void addEdge(const Edge* e);
	void removeEdge(const Edge* e);

	void printEdgesTo(const size_t genome) const;
	
	const Edge* getBestEdge(const size_t genome) const;
// 	Edge* getEdgeTo(Anchor* a) const;
// 	bool hasEdgeTo(Anchor* a) const;

	bool hasEdge(const Edge* e) const;

	bool hasEdgesTo(const size_t genome) const;

	// Filter the edges to GENOME by only keeping those that involve
	// anchors in between START and END
	void filterEdges(const size_t genome,
					 const Anchor* start,
					 const Anchor* end);

	// Filter the edges to GENOME by only keeping those that involve
	// anchors in between START and END
	void filterEdges(const size_t genome,
					 const Anchor* start1, const Anchor* end1,
					 const Anchor* start2, const Anchor* end2);

	void filterEdges(const vector< pair<Anchor*, Anchor*> >& ints1,
					 const vector< pair<Anchor*, Anchor*> >& ints2,
					 const Mask& filterMask);

	void filterEdges(const vector< pair<Anchor*, Anchor*> >& ints,
					 const Mask& filterMask);

	void removeAllEdges();
	void removeNonCliqueEdges();

	bool isInClique() const;
	bool isInRun() const;

	Clique* getClique() const;
	Run* getRun() const;
	
	Anchor* nextAnchor(bool forward = true) const;
	Anchor* nextAnchorInClique(bool forward = true) const;
	Anchor* nextAnchorInRun(bool forward = true) const;

	Clique* nextClique(bool forward = true) const;
	Run* nextRun(bool forward = true) const;

	bool operator<(const Anchor& other) const;
	friend ostream& operator<<(ostream& strm, const Anchor& a);

	void writeAnchorLine(ostream& strm) const;
	
	// Flip the strand of this anchor (coordinates stay the same)
	void flip();

	// Unflip this anchor if it has previously been flipped
	void unflip();

	// Changes the coordinates and strand of this anchor as if the
	// chromosome it is on is flipped
	void reverse();

	// Shift the coordinates of this anchor by OFFSET
	void shift(const GenomicDist offset);

	// Remove all edges incident to anchors in between this and OTHER
	void removeEdgesBetween(const Anchor* other);

	// Remove all edges incident to anchors in between this and the
	// end of this anchor's chromosome indicated by FORWARD
	void removeEdgesToChromEnd(const bool forward);
	
	void createClique(Clique& c, const Mask& genomes);
	
private:

	string name;
	Chromosome* chrom;
	char strand;
	GenomicDist start;
	GenomicDist end;
	size_t isCoding;

	bool flipped;
	
	Anchor* next;
	Anchor* prev;

	Clique* clique;

	vector< vector<const Edge*> > edgeVects;

	bool marked;
};

struct AnchorSorter {
	bool operator()(const Anchor* a1, const Anchor* a2) const {
		return (*a1) < (*a2);
	}
};

// Implementation of short Anchor members

inline GenomicDist Anchor::getStart() const { return start; }
inline GenomicDist Anchor::getEnd() const { return end; }
inline const string& Anchor::getName() const { return name; }
inline Chromosome* Anchor::getChrom() const { return chrom; }
inline Genome* Anchor::getGenome() const { return chrom->getGenome(); }
inline char Anchor::getStrand() const { return strand; }

inline size_t Anchor::getGenomeNum() const {
	return chrom->getGenome()->getNum();
}

inline string Anchor::getUniqueName() const {
	return chrom->getUniqueName() + "_" + name;
}

inline GenomicDist Anchor::getDistanceToChromEnd(const bool forward) const {
	return (forward ? chrom->getLength() - end : start);
}

inline int Anchor::getNumEdges(const size_t genome) const {
	return edgeVects[genome].size();
}

inline Clique* Anchor::getClique() const { return clique; }

inline Run* Anchor::getRun() const {
	return clique != NULL ? clique->getRun() : NULL;
}

inline GenomicDist Anchor::getGenomeStart() const {
	return chrom->getGenomeStart() + start;
}

inline GenomicDist Anchor::getGenomeEnd() const {
	return chrom->getGenomeStart() + end;
}

inline GenomicDist Anchor::getGenomeMiddle() const {
	return chrom->getGenomeStart() + (start + end) / 2;
}

inline GenomicDist Anchor::getLength() const { return end - start; }

inline GenomicDist Anchor::getDistanceTo(Anchor* a, const bool forward) const {
	assert(chrom == a->getChrom());
	return (forward ? a->getStart() - end : start - a->getEnd());
}

inline bool Anchor::isForward() const { return strand == '+'; }
inline bool Anchor::isMarked() const { return marked; }
inline bool Anchor::isInClique() const { return clique != NULL; }
inline bool Anchor::isInRun() const {
	return clique != NULL && clique->isInRun();
}

inline void Anchor::setNext(Anchor* a) { next = a; }	
inline void Anchor::setPrev(Anchor* a) { prev = a; }
inline void Anchor::setClique(Clique* c) { clique = c; }
inline void Anchor::setChrom(Chromosome* c) { chrom = c; }
inline void Anchor::setMarked(const bool state) { marked = state; }

inline const Edge* Anchor::getBestEdge(size_t genome) const {
	return edgeVects[genome].front();
}

inline bool Anchor::hasEdgesTo(const size_t genome) const {
	return !edgeVects[genome].empty();
}

inline Anchor* Anchor::nextAnchor(const bool forward) const {
	return forward ? next : prev;
}
inline Clique* Anchor::nextClique(const bool forward) const {
	Anchor* nextInClique = nextAnchorInClique(forward);
	return nextInClique != NULL ? nextInClique->getClique() : NULL;
}
inline Run* Anchor::nextRun(const bool forward) const {
	Anchor* nextInRun = nextAnchorInRun(forward);
	return nextInRun != NULL ? nextInRun->getRun() : NULL;
}

inline void Anchor::flip() {
	strand = (strand == '+' ? '-' : '+');
	flipped = !flipped;
}

inline void Anchor::unflip() {
	if (flipped) { flip(); }
}

inline void Anchor::shift(const GenomicDist offset) {
	start += offset;
	end += offset;
}

inline bool Anchor::operator<(const Anchor& a) const {
	return getGenomeStart() < a.getGenomeStart();
}

#endif // __ANCHOR_HH__
