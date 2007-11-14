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

#include "clique.hh"
#include "mask.hh"
#include "anchor.hh"
#include "genome.hh"

Clique::Clique() :
	mask(0),
	anchors(Genome::getNumGenomes(), NULL),
	run(NULL),
	keep(false)
{}

Clique::~Clique() {
	removeAnchorPtrs();
}

void Clique::keepClique() { keep = true; }

bool Clique::isKept() const { return keep; }

void Clique::flip() {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (mask[g]) {
			anchors[g]->flip();
		}
	}
}

void Clique::setMask(const Mask& newMask) {
	mask &= newMask;
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (not mask[g] and anchors[g] != NULL) {
			assert(anchors[g]->isInClique());
			anchors[g]->setClique(NULL);
			anchors[g] = NULL;
		}
	}
}

void Clique::setRun(Run* r) {
	run = r;
}

void Clique::addAnchor(Anchor* a) {
	anchors[a->getGenome()->getNum()] = a;
	mask.set(a->getGenome()->getNum());
}

void Clique::removeAnchor(Anchor* a) {
	assert(anchors[a->getGenome()->getNum()] == a);
	anchors[a->getGenome()->getNum()] = NULL;
	mask.reset(a->getGenome()->getNum());
	if (a->getClique() == this) {
		a->setClique(NULL);
	}
}

void Clique::removeAnchor(size_t g) {
	assert(hasGenome(g));
	removeAnchor(anchors[g]);
}

void Clique::removeAnchorPtrs() {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (anchors[g] != NULL) {
			if (anchors[g]->getClique() == this) {
				anchors[g]->setClique(NULL);
			}
		}
	}
}

void Clique::useClique() {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (mask[g]) {
			if (anchors[g]->isInClique()) {
				assert(anchors[g]->getClique() != this);
				delete anchors[g]->getClique();
			}
			anchors[g]->setClique(this);
			anchors[g]->removeNonCliqueEdges();
		}
	}
}

Mask Clique::getMask() const {
	return mask;
}

Anchor* Clique::getAnchor(size_t genome) const {
	return anchors[genome];
}

Anchor* Clique::getAnchor(Genome* g) const {
	return anchors[g->getNum()];
}

Run* Clique::getRun() const {
	return run;
}

size_t Clique::getSize() const {
	return mask.count();
}

int Clique::getWeight(const Mask& weightMask) const {
	return 0;
}

Clique* Clique::nextClique(size_t g, bool forward) {
	return anchors[g]->nextClique(forward);
}

Run* Clique::nextRun(size_t g, bool forward) {
	return anchors[g]->nextRun(forward);
}

bool Clique::isFirstClique(const size_t g) const {
	return anchors[g]->nextClique(false) == NULL;
}

bool Clique::isLastClique(const size_t g) const {
	return anchors[g]->nextClique(true) == NULL;
}

bool Clique::isEndClique(size_t g) const {
	return isFirstClique(g) or isLastClique(g);
}

bool Clique::isInRun() const {
	return run != NULL;
}

bool Clique::hasGenome(size_t g) const {
	return mask[g];
}

bool Clique::hasGenomes(const Mask& genomeMask) const {
	return (genomeMask & mask) == genomeMask;
}

void Clique::removeCliquesInBetween(const Clique* other,
									const Mask removeMask) {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (removeMask[g]) {
			Anchor* currAnchor = this->getAnchor(g);
			Anchor* nextAnchor = other->getAnchor(g);
			
			// Continue if anchors are not on the same chromosome (due
			// to draft genome).
			if (currAnchor->getChrom() != nextAnchor->getChrom()) {
				continue;
			}
			
			bool forward = *currAnchor < *nextAnchor;
			
			while (1) {
				// Walk towards the next clique's anchor
				currAnchor = currAnchor->nextAnchorInClique(forward);
				
				assert(currAnchor != NULL);
				
				// Check if we've arrived at the next clique's anchor
				if (currAnchor == nextAnchor) {
					break;
				} // If not, remove the clique at this anchor
				else {
					assert(!currAnchor->isInRun());
					delete currAnchor->getClique();
				}
			}
		}
	}
}

void Clique::removeEdgesInBetween(const Clique* other,
								  const size_t genome) {
	assert(hasGenome(genome) && other->hasGenome(genome));
	getAnchor(genome)->removeEdgesBetween(other->getAnchor(genome));
}

ostream& operator<<(ostream& strm, const Clique& c) {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (c.hasGenome(g)) {
			strm << (*c.getAnchor(g)) << '\t';
		}
	}
	strm << endl;
	return strm;
}

void Clique::print() const {
	std::cerr << *this << '\n';
}

bool Clique::isConnected() const {
	size_t g = firstInMask(mask);
	Clique c;
	getAnchor(g)->createClique(c, mask);
	return c == *this;
}

bool operator==(const Clique& c1, const Clique& c2) {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (c1.getAnchor(g) != c2.getAnchor(g)) {
			return false;
		}
	}
	return true;
}
