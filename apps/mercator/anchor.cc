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

#include "anchor.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "edge.hh"
#include "clique.hh"
#include "run.hh"

Anchor::Anchor(const string& name,
			   Chromosome* chrom,
			   const char strand,
			   const GenomicDist start,
			   const GenomicDist end,
			   const size_t isCoding) :
	name(name),
	chrom(chrom),
	strand(strand),
	start(start),
	end(end),
	isCoding(isCoding),

	flipped(false),

	next(NULL),
	prev(NULL),

	clique(NULL),

	edgeVects(Genome::getNumGenomes()),

	marked(false)
{}

void Anchor::writeAnchorLine(ostream& strm) const {
	strm << name << '\t'
		 << chrom->getName() << '\t'
		 << (flipped ? (strand == '+' ? '-' : '+') : strand) << '\t'
		 << start << '\t'
		 << end << '\t'
		 << isCoding << '\n';
}

void Anchor::print() const {
	std::cerr << name << ' ' << *this << '\n';
}

void Anchor::printEdgesTo(const size_t genome) const {
	for (size_t i = 0; i < edgeVects[genome].size(); ++i) {
		std::cerr << *(edgeVects[genome][i]) << '\n';
	}
}

size_t Anchor::getNameNum() const { return atoi(name.c_str()); }

bool Anchor::isRepetitive(size_t repeatNum, float repeatPct) const {
	assert(repeatNum > 1 && repeatPct <= 1.0 && repeatPct >= 0.0);

	// Check for repetitiveness in each genome
	vector< vector<const Edge*> >::const_iterator it;
	for (it = edgeVects.begin(); it != edgeVects.end(); ++it) {
		const vector<const Edge*>& edges = *it;
		// Definitely not repetitive if there is only one edge
		if (edges.size() <= 1) {
			continue;
		} // If the top two edges have the same score then it is repetitive
		else if (edges[0]->getScore() == edges[1]->getScore()) {
			return true;
		} // If the number of edges is less than REPEATNUM we can not
		  // be repetitive
		else if (edges.size() < repeatNum) {
			continue;
		} // Check the number of edges that have scores >= REPEATPCT * MAXSCORE
	    else {
			Edge::ScoreType rptScore = static_cast<Edge::ScoreType>(edges.front()->getScore() * repeatPct);
			size_t rptNum = 0;
			vector<const Edge*>::const_iterator eit;
			for (eit = edges.begin();
				 eit != edges.end() && (*eit)->getScore() >= rptScore;
				 ++eit, ++rptNum);

			if (rptNum >= repeatNum) {
				return true;
			} else {
				continue;
			}
		}
	}

	// If not found to be repetitive with respect to any genome, return false
	return false;
}

int Anchor::getMaxEdges() const {
	int max = 0;
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (max < getNumEdges(g)) {
			max = getNumEdges(g);
		}
	}
	return max;
}

void Anchor::removeEdge(const Edge* e) {
	if (isInClique() && clique->hasGenome(e->getOtherGenomeNum(this))) {
		std::cerr << "Found problem\n";
	}
	
	size_t genome = e->getOtherGenomeNum(this);

	// Find the edge in the edge vector
	std::vector<const Edge*>::iterator pos = 
		find(edgeVects[genome].begin(), edgeVects[genome].end(), e);

	// Assert that this edge actually exists
	assert(pos != edgeVects[genome].end());

	// Remove the edge from the edge vector
	edgeVects[genome].erase(pos);
}

ostream& operator<<(ostream& strm, const Anchor& a) {
	return strm << a.getChrom()->getName()
				<< a.getStrand()
				<< ':'
				<< a.getStart()
				<< '-'
				<< a.getEnd();
}

void Anchor::filterEdges(const size_t genome,
						 const Anchor* start,
						 const Anchor* end) {
// 	std::cerr << "Filtering edges between: ";
// 	if (start) {
// 		std::cerr << *start;
// 	} else {
// 		std::cerr << "-";
// 	}
// 	std::cerr << " and ";
// 	if (end) {
// 		std::cerr << *end;
// 	} else {
// 		std::cerr << "-";
// 	}
// 	std::cerr << " for anchor: " << *this << "\n";
	
	assert((start != NULL && genome == start->getGenome()->getNum()) ||
		   (end != NULL && genome == end->getGenome()->getNum()));
	// Step through edges and remove those that are not incident to
	// anchors in between START and END in GENOME
	vector<const Edge*>::iterator curr;
	vector<const Edge*>::iterator newend;
	for (curr = edgeVects[genome].begin(), newend = curr;
		 curr != edgeVects[genome].end(); ++curr) {
		Anchor* other = (*curr)->getOtherAnchor(this);
		// Check if this edge goes to an anchor between START and END
		// If START or END is NULL, we ignore that end point
		if ((start == NULL || (*start) < (*other)) &&
			(end == NULL || (*other) < (*end))) {
			(*newend) = (*curr); // push forward past removed edges
			++newend;
		} else {
// 			std::cerr << "Removing edge: " << **curr << "\n";
			other->removeEdge(*curr); // Remove edge from the other anchor
		}
	}

	// Truncate the edge vector to only include the edges that passed
	// the filter
	edgeVects[genome].erase(newend, edgeVects[genome].end());
}

void Anchor::filterEdges(const size_t genome,
						 const Anchor* start1, const Anchor* end1,
						 const Anchor* start2, const Anchor* end2) {
	// Step through edges and remove those that are not incident to
	// anchors in between START1 and END1 or between START2 and END2
	// in GENOME
	vector<const Edge*>::iterator curr;
	vector<const Edge*>::iterator newend;
	for (curr = edgeVects[genome].begin(), newend = curr;
		 curr != edgeVects[genome].end(); ++curr) {
		Anchor* other = (*curr)->getOtherAnchor(this);
		// Check if this edge goes to an anchor between START and END
		// If START or END is NULL, we ignore that end point
		if (((start1 == NULL || (*start1) < (*other))
			 && (end1 == NULL || (*other) < (*end1))
			 && (start1 != end1)) ||
			((start2 == NULL || (*start2) < (*other))
			 && (end2 == NULL || (*other) < (*end2))
			 && (start2 != end2))) {
			(*newend) = (*curr); // push forward past removed edges
			++newend;
		} else {
			// 			std::cerr << "Removing edge: " << **curr << "\n";
			other->removeEdge(*curr); // Remove edge from the other anchor
		}
	}
	
	// Truncate the edge vector to only include the edges that passed
	// the filter
	edgeVects[genome].erase(newend, edgeVects[genome].end());
}

void Anchor::filterEdges(const vector< pair<Anchor*, Anchor*> >& ints1,
						 const vector< pair<Anchor*, Anchor*> >& ints2,
						 const Mask& filterMask) {

	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (filterMask[g]) {
			this->filterEdges(g,
							  ints1[g].first, ints1[g].second,
							  ints2[g].first, ints2[g].second);
		}
	}
}

void Anchor::filterEdges(const vector< pair<Anchor*, Anchor*> >& ints,
						 const Mask& filterMask) {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (filterMask[g]) {
			this->filterEdges(g, ints[g].first, ints[g].second);
		}
	}
}

void Anchor::removeAllEdges() {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		vector<const Edge*>::iterator ePos;
		for (ePos = edgeVects[g].begin(); ePos != edgeVects[g].end(); ++ePos) {
			(*ePos)->getOtherAnchor(this)->removeEdge(*ePos);
		}
		edgeVects[g].clear();
	}
}

void Anchor::removeNonCliqueEdges() {
	assert(isInClique());
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (not clique->hasGenome(g) or g == getGenomeNum()) {
			continue;
		}
		const Edge* cliqueEdge = NULL;
		vector<const Edge*>::iterator ePos;
		for (ePos = edgeVects[g].begin(); ePos != edgeVects[g].end(); ++ePos) {
			Anchor* other = (*ePos)->getOtherAnchor(this);
			if (other == clique->getAnchor(g)) {
				cliqueEdge = *ePos;
			} else {
				other->removeEdge(*ePos);
			}
		}
		edgeVects[g].clear();

		if (cliqueEdge != NULL) {
			addEdge(cliqueEdge);
		}
	}
}

Anchor* Anchor::nextAnchorInClique(const bool forward) const {
	for (Anchor* curr = this->nextAnchor(forward);
		 curr != NULL;
		 curr = curr->nextAnchor(forward)) {
		if (curr->isInClique()) {
			return curr;
		}
	}
	return NULL;
}

Anchor* Anchor::nextAnchorInRun(const bool forward) const {
	for (Anchor* curr = this->nextAnchorInClique(forward);
		 curr != NULL;
		 curr = curr->nextAnchorInClique(forward)) {
		if (curr->isInRun()) {
			return curr;
		}
	}
	return NULL;
}

void Anchor::reverse() {
	// Flip the strand
	strand = (strand == '+' ? '-' : '+');

	// Flip the coordinates
	GenomicDist newStart = chrom->getLength() - end;
	end = chrom->getLength() - start;
	start = newStart;

	// Flip the next and prev pointers
	Anchor* temp = next;
	next = prev;
	prev = temp;
}

void Anchor::removeEdgesBetween(const Anchor* other) {
	assert(getChrom() == other->getChrom());

	bool forward = *this < *other;
	
	for (Anchor* curr = nextAnchor(forward);
		 curr != NULL && curr != other;
		 curr = curr->nextAnchor(forward)) {
		assert(!curr->isInClique() && !curr->isInRun());
		curr->removeAllEdges();
	}
}

void Anchor::removeEdgesToChromEnd(const bool forward) {
	for (Anchor* curr = nextAnchor(forward);
		 curr != NULL;
		 curr = curr->nextAnchor(forward)) {
		assert(!curr->isInClique() && !curr->isInRun());
		curr->removeAllEdges();
	}
}

bool Anchor::hasEdge(const Edge* e) const {
	size_t g = e->getOtherGenomeNum(this);
	return find(edgeVects[g].begin(), edgeVects[g].end(), e)
		!= edgeVects[g].end();
}

void Anchor::addEdge(const Edge* e) {
	size_t genome = e->getOtherGenomeNum(this);
	edgeVects[genome].push_back(e);
}


void Anchor::createClique(Clique& c, const Mask& genomes) {
	if (c.hasGenome(getGenomeNum())) {
		if (c.getAnchor(getGenomeNum()) != this) {
			std::cerr << "Clique does not include anchor!";
			assert(false);
		}
		return;
	}
	
	c.addAnchor(this);
	
	for (size_t i = 0; i < Genome::getNumGenomes(); ++i) {
		if (genomes[i] && i != getGenomeNum()) {
			if (hasEdgesTo(i)) {
				Anchor* a = getBestEdge(i)->getOtherAnchor(this);
				a->createClique(c, genomes);
			}
		}
	}
	
}
