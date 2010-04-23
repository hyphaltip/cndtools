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

#include "run.hh"

#include "util/stl.hh"
#include "anchor.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "clique.hh"
#include "mask.hh"

Run::Run() :
	mask(0),
	cliques(),
	starts(Genome::getNumGenomes()),
	ends(Genome::getNumGenomes()),
	visited(false),
	joined(false)
{}

Run::~Run() {
	unclaimCliques();
}

size_t Run::minRunLength = 0;

void Run::setMinRunLength(size_t num) { minRunLength = num; }

size_t Run::getMinRunLength() { return minRunLength; }

bool Run::isSignificant() const {
	return (cliques.size() >= minRunLength or
			std::find_if(cliques.begin(), cliques.end(),
						 std::mem_fun(&Clique::isKept)) != cliques.end());
}

Mask Run::getMask() const { return mask; }

const vector<Clique*>& Run::getCliques() const { return cliques; }

bool Run::wasVisited() const {
	return visited;
}

void Run::setVisited(bool visited) {
	this->visited = visited;
}

bool Run::wasJoined() const {
	return joined;
}

ostream& operator<<(ostream& strm, const Run& r) {
	for (size_t i = 0; i < r.cliques.size(); ++i) {
		strm << *(r.cliques[i]);
	}
	return strm;
}

bool Run::hasAnchor(size_t num, size_t genome) {
	std::vector<Clique*>::iterator it;
	for (it = cliques.begin(); it != cliques.end(); ++it) {
		if ((*it)->hasGenome(genome) &&
			(*it)->getAnchor(genome)->getNameNum() == num) {
			return true;
		}
	}
	return false;
}

void Run::print() const {
	std::cerr << *this;
}

bool Run::isAssembled(const size_t genome) const {
	Chromosome* chrom = NULL;
	for (size_t i = 0; i < cliques.size(); ++i) {
		if (cliques[i]->hasGenome(genome)) {
			if (chrom == NULL) {
				chrom = cliques[i]->getAnchor(genome)->getChrom();
			} else if (chrom != cliques[i]->getAnchor(genome)->getChrom()) {
				return true;
			}
		}
	}
	// If we made it out of the loop, all chromosomes were the same
	return false;
}

bool Run::isLeftChromEnd(const size_t genome) const {
	if (isLeftForward(genome)) {
		return getLeftClique(genome)->isFirstClique(genome);
	} else {
		return getLeftClique(genome)->isLastClique(genome);
	}
}

bool Run::isRightChromEnd(const size_t genome) const {
	if (isRightForward(genome)) {
		return getRightClique(genome)->isLastClique(genome);
	} else {
		return getRightClique(genome)->isFirstClique(genome);
	}
}

bool Run::isLeftOf(const Run& other,
				   const size_t minAdjacent,
				   const GenomicDist maxDist,
				   const bool strict) const {
	size_t adjacent = 0;
	
	// Find genomes shared between the runs
	Mask shared = mask & other.getMask();

	// If the two runs do not share at least minAdjacent genomes,
	// return false
	if (shared.count() < minAdjacent) {
		return false;
	}

	// Check each shared genome to make sure the ends are consistent
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (not shared[g]) {
			continue;
		}

		Anchor* right = getRightAnchor(g);
		Anchor* left = other.getLeftAnchor(g);

		if (Genome::getGenome(g)->isDraft() and
			(right->getChrom() != left->getChrom())) {

			if ((right->getDistanceToChromEnd(isRightForward(g)) +
				 left->getDistanceToChromEnd(not other.isLeftForward(g)))
				> maxDist) {
				//cerr << "isLeftOf: No join because distances to chrom ends "
				//	 << "add up to more than maxdist" << endl;
				return false;
			}
			
			if (strict and
				(not isRightChromEnd(g) or not other.isLeftChromEnd(g))) {
				//cerr << "isLeftOf: No join because cliques aren't at chrom ends"
				//	 << endl;
				return false;
			}
			
		} else {
			if (right->getChrom() != left->getChrom()) {
				//cerr << "isLeftOf: chroms don't match in genome " << g << endl;
				return false;
			}

			// The runs are adjacent
			++adjacent;

			// This code can be reduced!
			if (other.isLeftForward(g)) {
				if (not isRightForward(g)) {
					//cerr << "isLeftOf: orientations don't match in genome "
					//	 << g << endl;
					return false;
				}
				if (*left < *right) {
					//cerr << "isLeftOf: order of anchors does not match "
					//	 << "orientation in genome " << g << endl;
					return false;
				}
				if (right->getDistanceTo(left, true) > maxDist) {
					//cerr << "isLeftOf: distance between anchors is too great "
					//	 << "in genome " << g << endl;
					return false;
				}
			} else {
				if (isRightForward(g)) {
					//cerr << "isLeftOf: orientations don't match in genome "
					//	 << g << endl;
					return false;
				}
				if (*right < *left) {
					//cerr << "isLeftOf: order of anchors does not match "
					//	 << "orientation in genome " << g << endl;
					return false;
				}
				if (right->getDistanceTo(left, false) > maxDist) {
					//cerr << "isLeftOf: distance between anchors is too great "
					//	 << "in genome " << g << endl;
					return false;
				}
			}
		}
	}

	if (adjacent < minAdjacent) {
		//cerr << "isLeftOf: too few adjacencies" << endl;
		return false;
	} else {
		return true;
	}
}

Run* Run::nextRun(const size_t genome,
				  const bool forward) const {
	return (isForward(genome) == forward ?
			ends[genome]->nextRun(genome, forward) :
			starts[genome]->nextRun(genome, forward));
}

Run* Run::leftRun(size_t genome) const {
	return starts[genome]->nextRun(genome, !isLeftForward(genome));
}

Run* Run::rightRun(size_t genome) const {
	return ends[genome]->nextRun(genome, isRightForward(genome));
}

Clique* Run::getLastClique() const {
	return cliques.back();
}

Clique* Run::getFirstClique() const {
	return cliques.front();
}

bool Run::isEmpty() const {
	return cliques.empty();
}

void Run::addCliqueToRight(Clique* c) {
	assert(c != NULL);
		
	// Add to vector of cliques
	cliques.push_back(c);

	// Update run mask to include genomes in this clique
	mask |= c->getMask();

	assert(starts.size() != 0);
	// Update end (and possibly start) vectors with this clique
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (c->hasGenome(g)) {
			if (starts[g] == NULL) {
				starts[g] = c;
			}
			ends[g] = c;
		}
	}
}

void Run::addRunToRight(Run* r) {
	vector<Clique*>::iterator c;
	for (c = r->cliques.begin(); c != r->cliques.end(); ++c) {
		addCliqueToRight(*c);
	}
}

void Run::split() {
	unclaimCliques();
	vector<Clique*>::iterator it;
	for (it = cliques.begin(); it != cliques.end(); ++it) {
		Run* r = new Run();
		r->addCliqueToRight(*it);
		r->claimCliques();
	}
	reset();
}

Run* Run::splitRunRight(const size_t genome) {
	// Find the index of the rightmost clique that is not on the right
	// chromosome for GENOME
	size_t splitIndex = cliques.size() - 1;
	while (!cliques[splitIndex]->hasGenome(genome) ||
		   cliques[splitIndex]->getAnchor(genome)->getChrom() ==
		   getRightChrom(genome)) {
		--splitIndex;
		assert(splitIndex >= 0);
	}

	// Make a copy of the cliques originally in this run
	vector<Clique*> cliquesToSplit(cliques.begin(), cliques.end());

	// Clear out the cliques in this run
	unclaimCliques();
	reset();

	// Make a new run for the cliques to split off to the right
	Run* rightRun = new Run();

	// Copy the cliques into the appropriate runs
	for (size_t i = 0; i < cliquesToSplit.size(); ++i) {
		if (i <= splitIndex) {
			addCliqueToRight(cliquesToSplit[i]);
		} else {
			rightRun->addCliqueToRight(cliquesToSplit[i]);
		}
	}

	// Make the two split runs claim their cliques
	claimCliques();
	rightRun->claimCliques();

	return rightRun;
}

int Run::getWeight() const {
	int weight = 0;
	for (size_t i = 0; i < cliques.size(); ++i) {
		weight += cliques[i]->getWeight(mask);
	}
	return weight;
}
	
void Run::reset() {
	mask.reset();
	cliques.clear();
	fill(starts.begin(), starts.end(), static_cast<Clique*>(NULL));
	fill(ends.begin(), ends.end(), static_cast<Clique*>(NULL));
}

void Run::flip() {
	std::reverse(cliques.begin(), cliques.end());
	starts.swap(ends);
	for_each(cliques.begin(), cliques.end(), mem_fun(&Clique::flip));
}
	
bool Run::hasGenomes(const Mask& genomeMask) const {
	return (genomeMask & mask) == genomeMask;
}

bool Run::hasGenome(const size_t genome) const {
	return mask[genome];
}
	
size_t Run::getLength() const {
	return cliques.size();
}

bool Run::isForward(const size_t genome) const {
	assert(hasGenome(genome));
	assert(!isAssembled(genome));
	return starts[genome]->getAnchor(genome)->isForward();
}

bool Run::isLeftForward(const size_t genome) const {
	return starts[genome]->getAnchor(genome)->isForward();
}
	
bool Run::isRightForward(const size_t genome) const {
	return ends[genome]->getAnchor(genome)->isForward();
}

pair<Anchor*, Anchor*>
Run::getLeftInterRunAnchorInterval(const size_t genome) const {
	if (!hasGenome(genome)) {
		return make_pair<Anchor*, Anchor*>(NULL, NULL);
	}
	Anchor* leftAnchor = this->getLeftAnchor(genome);
	Anchor* nextAnchor = leftAnchor->nextAnchorInRun(!leftAnchor->isForward());
	return (leftAnchor->isForward() ?
			make_pair(nextAnchor, leftAnchor) :
			make_pair(leftAnchor, nextAnchor));
}

pair<Anchor*, Anchor*>
Run::getRightInterRunAnchorInterval(const size_t genome) const {
	if (!hasGenome(genome)) {
		return make_pair<Anchor*, Anchor*>(NULL, NULL);
	}
	Anchor* rightAnchor = this->getRightAnchor(genome);
	Anchor* nextAnchor = rightAnchor->nextAnchorInRun(rightAnchor->isForward());
	return (rightAnchor->isForward() ?
			make_pair(rightAnchor, nextAnchor) :
			make_pair(nextAnchor, rightAnchor));
}

void
Run::getLeftInterRunAnchorIntervals(vector< pair<Anchor*, Anchor*> >& ints) const {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		ints[g] = this->getLeftInterRunAnchorInterval(g);
	}
}

void
Run::getRightInterRunAnchorIntervals(vector< pair<Anchor*, Anchor*> >& ints) const {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		ints[g] = this->getRightInterRunAnchorInterval(g);
	}
}

Anchor* Run::getLeftAnchor(const size_t genome) const {
	return getLeftClique(genome)->getAnchor(genome);
}

Anchor* Run::getRightAnchor(const size_t genome) const {
	return getRightClique(genome)->getAnchor(genome);
}

Clique* Run::getLeftClique(const size_t genome) const {
	assert(hasGenome(genome));
	return starts[genome];
}

Clique* Run::getRightClique(const size_t genome) const {
	assert(hasGenome(genome));
	return ends[genome];
}

void Run::getCliques(vector<Clique*>& gcliques, size_t g) const {
	util::stl::copy_if(cliques.begin(), cliques.end(),
					 std::back_inserter(gcliques),
					 std::bind2nd(std::mem_fun(&Clique::hasGenome), g));
}

Chromosome* Run::getLeftChrom(const size_t genome) const {
	return getLeftAnchor(genome)->getChrom();
}

Chromosome* Run::getRightChrom(const size_t genome) const {
	return getRightAnchor(genome)->getChrom();
}

Chromosome* Run::getChrom(const size_t genome) const {
	assert(hasGenome(genome));
	assert(!isAssembled(genome));
	return getLeftChrom(genome);
}

size_t Run::getNumCliques(const size_t genome) const {
	assert(hasGenome(genome));	
	size_t numCliques = 0;
	for (Anchor* a = getLeftAnchor(genome);
		 a != getRightAnchor(genome);
		 a = a->nextAnchor(isForward(genome))) {
		if (a->isInClique()) {
			++numCliques;
		}
	}
	return numCliques;
}

size_t Run::getNumAnchors(const size_t genome) const {
	assert(hasGenome(genome));
	size_t numAnchors = 1;
	for (Anchor* a = getLeftAnchor(genome);
		 a != getRightAnchor(genome);
		 a = a->nextAnchor(isForward(genome))) {
		++numAnchors;
	}
	return numAnchors;
}

GenomicDist Run::leftCoord(const size_t genome) const {
	assert(hasGenome(genome));
	return (isLeftForward(genome) ?
			getLeftAnchor(genome)->getStart() :
			getLeftAnchor(genome)->getEnd());
}

GenomicDist Run::rightCoord(const size_t genome) const {
	assert(hasGenome(genome));
	return (isRightForward(genome) ?
			getRightAnchor(genome)->getEnd() :
			getRightAnchor(genome)->getStart());
}

GenomicDist Run::leftSplitCoord(const size_t genome) const {
	assert(hasGenome(genome));
	Run* lr = leftRun(genome);
	assert(lr != this);
	if (lr != NULL && lr->getChrom(genome) == this->getChrom(genome)) {
		if (isForward(genome)) {
			return (lr->endCoord(genome) + startCoord(genome)) / 2;
		} else {
			return (lr->startCoord(genome) + endCoord(genome)) / 2;
		}
	} else {
		if (isForward(genome)) {
			return max(this->leftCoord(genome) - 1000, GenomicDist(0));
		} else {
			return min(this->leftCoord(genome) + 1000,
					   getChrom(genome)->getLength());
		}
	}
}

GenomicDist Run::rightSplitCoord(const size_t genome) const {
	assert(hasGenome(genome));
	Run* rr = rightRun(genome);
	assert(rr != this);
	if (rr != NULL && rr->getChrom(genome) == this->getChrom(genome)) {
		if (isForward(genome)) {
			return (rr->startCoord(genome) + endCoord(genome)) / 2;
		} else {
			return (rr->endCoord(genome) + startCoord(genome)) / 2;
		}
	} else {
		if (isForward(genome)) {
			return min(this->rightCoord(genome) + 1000,
					   getChrom(genome)->getLength());
		} else {
			return max(this->rightCoord(genome) - 1000, GenomicDist(0));
		}
	}
}

GenomicDist Run::startCoord(const size_t genome) const {
	assert(hasGenome(genome));
	return min(getLeftAnchor(genome)->getStart(),
			   getRightAnchor(genome)->getStart());
}

GenomicDist Run::endCoord(const size_t genome) const {
	assert(hasGenome(genome));
	return max(getLeftAnchor(genome)->getEnd(),
			   getRightAnchor(genome)->getEnd());
}

GenomicDist Run::start(const size_t genome) const {
	assert(hasGenome(genome));
	return min(getLeftAnchor(genome)->getGenomeStart(),
			   getRightAnchor(genome)->getGenomeStart());
}

GenomicDist Run::end(const size_t genome) const {
	assert(hasGenome(genome));
	return max(getLeftAnchor(genome)->getGenomeEnd(),
			   getRightAnchor(genome)->getGenomeEnd());
}

GenomicDist Run::getCoverage(const size_t genome) const {
	assert(hasGenome(genome));
	return endCoord(genome) - startCoord(genome);
}

void Run::unclaimCliques() {
	vector<Clique*>::iterator c;
	for (c = cliques.begin(); c != cliques.end(); ++c) {
		if ((*c)->getRun() == this) {
			(*c)->setRun(NULL);
		}
	}
}

void Run::claimCliques() {
	vector<Clique*>::iterator c;
	for (c = cliques.begin(); c != cliques.end(); ++c) {
		if ((*c)->isInRun()) {
			delete (*c)->getRun();
		}
		(*c)->setRun(this);
	}
}

void Run::deleteCliques() {
	for (vector<Clique*>::iterator c = cliques.begin();
		 c != cliques.end(); ++c) {
		delete (*c);
	}
}	

void Run::dump(ostream& strm) const {
	vector<Clique*>::const_iterator c;
	for (c = cliques.begin(); c != cliques.end(); ++c) {
		assert((*c)->getSize() > 0);
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if (g > 0) {
				strm << '\t';
			}

			if ((*c)->hasGenome(g)) {
				strm << (*c)->getAnchor(g)->getName();
			} else {
				strm << "NA";
			}
		}
		strm << endl;
	}
}

void Run::printMapLine(ostream& strm,
					   bool extend) const {
	strm << getNum();

	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		strm << '\t';
		if (hasGenome(g)) {
			GenomicDist start = (extend ? leftSplitCoord(g) : leftCoord(g));
			GenomicDist end = (extend ? rightSplitCoord(g) : rightCoord(g));
			if (start > end) {
				std::swap(start, end);
			}
			strm << getChrom(g)->getName()
				 << '\t'
				 << start
				 << '\t'
				 << end
				 << '\t'
				 << (isForward(g) ? '+' : '-');
		} else {
			strm << "NA\tNA\tNA\tNA";
		}
	}

	strm << endl;
}

void Run::printPairwiseHits(ostream& strm) const {
	vector<Clique*>::const_iterator c;
	for (c = cliques.begin(); c != cliques.end(); ++c) {
		assert((*c)->getSize() > 0);
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if ((*c)->hasGenome(g)) {
				for (size_t h = g + 1; h < Genome::getNumGenomes(); ++h) {
					if ((*c)->hasGenome(h)
						&& (*c)->getAnchor(g)->hasEdgesTo(h)) {
						strm << getNum() << '\t'
							 << Genome::getGenome(g)->getName() << '\t'
							 << (*c)->getAnchor(g)->getName() << '\t'
							 << Genome::getGenome(h)->getName() << '\t'
							 << (*c)->getAnchor(h)->getName() << endl;
					}
				}
			}
		}
	}
}

void Run::resetJoinFlags() {
	visited = false;
	joined = false;
	
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (mask[g]) {
			Run* next = leftRun(g);
			if (next != NULL && next->wasVisited()) {
				next->resetJoinFlags();
			}
		}
	}
}

Run* Run::joinTo(Run& other, const GenomicDist maxDist) {
	Run* joinedRun = new Run();
	bool joinPossible = canJoinTo(other, *joinedRun, maxDist);
	resetJoinFlags();

	if (joinPossible) {
		return joinedRun;
	} else {
		delete joinedRun;
		return NULL;
	}
}

bool Run::canJoinTo(Run& other, Run& joinedRun, const GenomicDist maxDist) {
	visited = true;

	if (this == &other) {
		joined = true;
		joinedRun.addRunToRight(this);
		//		cerr << "Target run hit, returned true" << endl;
		return true;
	} else if (!other.isLeftOf(*this)) {
		return false;
	}

	// Find genomes shared between runs
	Mask shared = mask & other.getMask();

	// Mask for indicating whether there was a successful traversal
	// towards the target run along each genome
	Mask success = 0;

	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (mask[g]) {

			// We are searching from right to left towards the target run
			Run* next = leftRun(g);

			// cerr << "Left run:" << endl << *next << endl;

			// If the run to the left has not already been visited,
			// make sure it is in the right orientation
			if (next != NULL && !next->wasVisited()) {
				// Reverse the right run if this is a draft genome and
				// the ends are not from the same chromosome or if the
				// orientations of the ends aren't the same
				if ((Genome::getGenome(g)->isDraft() &&
					 getLeftChrom(g) != next->getRightChrom(g)) ||
					isLeftForward(g) != next->isRightForward(g)) {
					next->flip();

					// Assert that the chroms and orientations are the
					// same after this check
					assert((getLeftChrom(g) == next->getRightChrom(g)) &&
						   (isLeftForward(g) == next->isRightForward(g)));
				}
			}

			// Check that the next run has the following properties:
			// 1 - Not NULL if this is not a draft genome
			// 2 - If it has been visited already, it has been joined
			// 3 - It is to the left of this run
			// 4 - It can be joined to the target run
			if ((next == NULL && Genome::getGenome(g)->isDraft()) ||
				(next != NULL &&
				 next->isLeftOf(*this, 1, maxDist) &&
				 (next->wasJoined() ||
				  (!next->wasVisited() &&
				   next->canJoinTo(other, joinedRun, maxDist))))) {
				success[g] = true;
			}
			// If this is a shared genome with the target run,
			// then this run can not be successfully joined, so
			// break out of the loop.
			else if (shared[g]) {
				break;
			} else {
				continue;
			}
		}
	}

	// Update the shared mask to be those genomes that are shared
	// between this run and the current joined run
	shared = mask & joinedRun.mask;
	
	// This run must have at least one shared genome and must have
	// success on all shared genomes
	if (shared.none() || (success & shared) != shared) {
		// cerr << "Can not join: no success: " << success << endl;
		return false;
	}

	// This run must also be able to strictly join the joined run at
	// the right of the joined run
	if (!joinedRun.isLeftOf(*this, 1, maxDist, true)) {
		// cerr << "Can not join: previously joined is not left of this" << endl;
		return false;
	}

	// If all is well by this point, we are ready to join this run to
	// the joined run, mark this run as joined, and return true
	// indicating success
	joinedRun.addRunToRight(this);
	joined = true;
	return true;
}

void Run::filterInterRunEdges() {
	assert(!this->wasVisited());
	vector< pair<Anchor*, Anchor*> > leftIntervals(Genome::getNumGenomes());
	vector< pair<Anchor*, Anchor*> > rightIntervals(Genome::getNumGenomes());

	// Filter once to the right and once to the left
	for (int flip = 0; flip <= 1; ++flip) {
		if (flip == 1) {
			this->flip();
		}
		
		// Filter anchors to the left
		this->getLeftInterRunAnchorIntervals(leftIntervals);
		
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if (!hasGenome(g)) {
				continue;
			}
			
			Run* left = leftRun(g);
			if (left == NULL) {
				Mask filterMask = this->getMask();
				filterMask.reset(g);
				for (Anchor* a = getLeftAnchor(g)->nextAnchor(!isLeftForward(g));
					 a != NULL;
					 a = a->nextAnchor(!isLeftForward(g))) {
					a->filterEdges(leftIntervals, filterMask);
				}
			} else if (!left->wasVisited()) {
				if ((left->getRightChrom(g) != this->getLeftChrom(g)) ||
					(left->isRightForward(g) != this->isLeftForward(g))) {
					left->flip();
					assert((left->getRightChrom(g) == this->getLeftChrom(g)) &&
						   (left->isRightForward(g) == this->isLeftForward(g)));
				}
				
				left->getRightInterRunAnchorIntervals(rightIntervals);
				
				Mask filterMask = this->getMask() | left->getMask();
				filterMask.reset(g);
				for (Anchor* a = leftIntervals[g].first->nextAnchor();
					 a != leftIntervals[g].second;
					 a = a->nextAnchor()) {
					a->filterEdges(leftIntervals, rightIntervals, filterMask);
				}
			}
			
		}
		
	}

	this->setVisited(true);
}

void Run::filterIntraRunEdges() {
	typedef vector<Clique*>::iterator CliqueIter;

	// First establish some outer limits
	vector<Anchor*> outerBegins(Genome::getNumGenomes(), NULL);
	vector<Anchor*> outerEnds(Genome::getNumGenomes(), NULL);

	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (!hasGenome(g)) {
			continue;
		}
		outerBegins[g] =
			getLeftAnchor(g)->nextAnchorInRun(!getLeftAnchor(g)->isForward());
		outerEnds[g] =
			getRightAnchor(g)->nextAnchorInRun(getRightAnchor(g)->isForward());
	}

	// Filter internal edges in each genome in this run
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (!hasGenome(g)) {
			continue;
		}
		
		// Find the first clique that includes G
		CliqueIter begin =
			std::find_if(cliques.begin(), cliques.end(),
						 std::bind2nd(std::mem_fun(&Clique::hasGenome), g));
		CliqueIter end = begin;
		
		assert(begin != cliques.end());
		
		while (true) {
			begin = end;
			// Find the next clique that includes G
			end =
				std::find_if(begin + 1, cliques.end(),
							 std::bind2nd(std::mem_fun(&Clique::hasGenome), g));

			// If no more cliques include G, then we are done
			if (end == cliques.end()) {
				break;
			}

			// Vectors to mark the intervals in the other genomes
			// where edges should be allowed from the current genome G
			vector<Anchor*> begins(Genome::getNumGenomes(), NULL);
			vector<Anchor*> ends(Genome::getNumGenomes(), NULL);

			// Masks indicating the genomes for which we have found
			// interval end points
			Mask beginMask((*begin)->getMask());
			Mask endMask((*end)->getMask());

			// Find the begin interval end points by stepping backward
			// in the clique vector and adding begin points from cliques
			// that have a path to the begin clique.  We keep track of
			// possible paths by updating the beginMask with genomes
			// that have cliques as part of the path.
			CliqueIter curr;
			for (curr = begin; curr >= cliques.begin(); --curr) {
				// Ignore this clique if does not share any genomes
				// with the current mask (meaning that it does not
				// have a path to the begin clique)
				if (((*curr)->getMask() & beginMask).any()) {
					// Update the begin points for any genomes that we
					// have not seen yet
					for (size_t f = 0; f < Genome::getNumGenomes(); ++f) {
						if (begins[f] == NULL && (*curr)->hasGenome(f)) {
							begins[f] = (*curr)->getAnchor(f);
						}
					}
					// Update the genomes that we have seen
					beginMask |= (*curr)->getMask();
					if (beginMask == mask) {
						break;
					}
				}
			}

			// Do the same as above, but for the end points.  Here we
			// step forward in the clique vector.
			for (curr = end; curr != cliques.end(); ++curr) {
				if (((*curr)->getMask() & endMask).any()) {
					for (size_t f = 0; f < Genome::getNumGenomes(); ++f) {
						if (ends[f] == NULL && (*curr)->hasGenome(f)) {
							ends[f] = (*curr)->getAnchor(f);
						}
					}
 					endMask |= (*curr)->getMask();
					if (endMask == mask) {
						break;
					}
				}
			}

			Anchor* beginAnchor = (*begin)->getAnchor(g);
			Anchor* endAnchor = (*end)->getAnchor(g);

			assert((beginAnchor->getChrom() != endAnchor->getChrom()) ||
				   (beginAnchor->isForward() == endAnchor->isForward()));

			// Remove edges incident to anchors in between begin and
			// end that are not incident to anchors within the
			// intervals of the other genomes discovered above
			for (size_t f = 0; f < Genome::getNumGenomes(); ++f) {
				if (f == g || !hasGenome(f)) {
					continue;
				}

				bool beginForward;
				bool endForward;
				
				// Add outer end points if needed
				if (begins[f] == NULL) {
					begins[f] = outerBegins[f];
					beginForward = getLeftAnchor(f)->isForward();
				} else {
					beginForward = begins[f]->isForward();
				}
					
				if (ends[f] == NULL) {
					ends[f] = outerEnds[f];
					endForward = getRightAnchor(f)->isForward();
				} else {
					endForward = ends[f]->isForward();
				}

				// If no limits can be established, do not filter
				if (begins[f] == NULL && ends[f] == NULL) {
					continue;
				}
				
				// If genome F is a draft genome, and the two end
				// points are on different contigs, do not filter any
				// edges.  NOTE: This is a behavior that we may want to
				// change at some point.
				if (Genome::getGenome(f)->isDraft() &&
					(begins[f] == NULL || ends[f] == NULL ||
					 begins[f]->getChrom() != ends[f]->getChrom())) {
					continue;
				}

				assert(beginForward == endForward);
				// Swap end points of interval in genome F if it is
				// not forward
				if (!beginForward) {
					// Check for cyclic run, and make sure it really is cyclic
					if (begins[f] != NULL && ends[f] != NULL &&
						*(begins[f]) < *(ends[f])) {
						assert(isAssembled(f) &&
							   (getLeftChrom(f) == getRightChrom(f)));
						continue;
					}
					std::swap(begins[f], ends[f]);
				} else {
					if (begins[f] != NULL && ends[f] != NULL &&
						*(ends[f]) < *(begins[f])) {
						assert(isAssembled(f) &&
							   (getLeftChrom(f) == getRightChrom(f)));
						continue;
					}
				}

				// Check if our anchors are on different contigs of a
				// draft genome
				if (beginAnchor->getChrom() != endAnchor->getChrom()) {
					// Filter to end of contig from begin anchor
					for (Anchor* a = beginAnchor->nextAnchor(beginAnchor->isForward());
						 a != NULL;
						 a = a->nextAnchor(beginAnchor->isForward())) {
						a->filterEdges(f, begins[f], ends[f]);
					}
					// Filter to end of contig from end anchor
					for (Anchor* a = endAnchor->nextAnchor(!endAnchor->isForward());
						 a != NULL;
						 a = a->nextAnchor(!endAnchor->isForward())) {
						a->filterEdges(f, begins[f], ends[f]);
					}
				} else {
					// Filter between begin and end anchors
					for (Anchor* a = beginAnchor->nextAnchor(beginAnchor->isForward());
						 a != endAnchor;
						 a = a->nextAnchor(beginAnchor->isForward())) {
						a->filterEdges(f, begins[f], ends[f]);
					}
				}
				
			}

		}
		
	}
}

void Run::removeInterveningEdges() {
	// Keep track of the last anchors visited in each genome
	vector<Anchor*> currAnchors(Genome::getNumGenomes(), NULL);

	// Step through the cliques from left to right and filter the
	// edges incident to the anchors in between the clique anchors
	vector<Clique*>::iterator it;
	for (it = cliques.begin(); it != cliques.end(); ++it) {
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if ((*it)->hasGenome(g)) {
				Anchor* curr = currAnchors[g];
				Anchor* next = (*it)->getAnchor(g);
				if (curr != NULL) {
					if (curr->getChrom() == next->getChrom()) {
						curr->removeEdgesBetween(next);
					} else {
						curr->removeEdgesToChromEnd(curr->isForward());
						next->removeEdgesToChromEnd(!next->isForward());
					}
				}
				currAnchors[g] = (*it)->getAnchor(g);
			}
		}
	}
}

bool RunLengthSorter::operator()(const Run* r1, const Run* r2) const {
	return r1->getLength() < r2->getLength();
}

bool RunSorter::operator()(const Run* r1, const Run* r2) const {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (r1->hasGenome(g)) {
			return !r2->hasGenome(g) || r1->start(g) < r2->start(g);
		} else if (r2->hasGenome(g)) {
			return false;
		}
	}
	assert(false); // We should always return before this point
	return false;
}

void orderRuns(vector<Run*>& runs) {
	vector<Run*>::iterator runPtr;
	for (runPtr = runs.begin(); runPtr != runs.end(); ++runPtr) {
		if (!(*runPtr)->isForward(firstInMask((*runPtr)->getMask()))) {
			(*runPtr)->flip();
		}
	}
	sort(runs.begin(), runs.end(), RunSorter());
}

void getRuns(vector<Run*>& runs) {
	// Make a mask with all bits set
	Mask included = 0;
	included.flip();

	getRuns(runs, included);
}

void getRuns(vector<Run*>& runs, const size_t g) {
	// Make a mask with just the bit for G set
	Mask included = 0;
	included[g] = true;
	getRuns(runs, included);
}

void getRuns(vector<Run*>& runs, const Mask& included) {
	set<Run*> runSet;

	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		if (!included[g]) {
			continue;
		}
		Genome* gen = Genome::getGenome(g);
		for (size_t c = 0; c < gen->getNumChroms(); ++c) {
			for (Clique* curr = gen->getChrom(c)->getFirstClique();
				 curr != NULL;
				 curr = curr->nextClique(g)) {
				if (curr->isInRun()) {
					runSet.insert(curr->getRun());
				}
			}
		}
	}

	copy(runSet.begin(), runSet.end(), back_inserter(runs));
}

void joinRuns(const size_t minAdjacent,
			  const GenomicDist maxDist) {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		Genome* gen = Genome::getGenome(g);

		// cerr << "Genome " << g << endl;

		for (size_t c = 0; c < gen->getNumChroms(); ++c) {
			Chromosome* chrom = gen->getChrom(c);

			// Start at the first (left-most) run in this chromosome
			Run* prevRun = chrom->getFirstRun();

			// If there are no runs in this chromosome, continue
			if (prevRun == NULL) {
				continue;
			}

			// Orient the run so that its right end is forward and on
			// this chromosome
			if (prevRun->getRightChrom(g) != chrom or
				not prevRun->isRightForward(g)) {
				prevRun->flip();
			}

			// If the right end is still not forward and on this
			// chromosome, then this must be an assembled run without
			// an end to extend off of on this chromosome.  Just
			// continue on.
			if (prevRun->getRightChrom(g) != chrom or
				not prevRun->isRightForward(g)) {
				assert(prevRun->isAssembled(g));
				continue;
			}

			while (true) {

				// Find the next run to the right of this one
				Run* nextRun = prevRun->rightRun(g);

				// Stop when we have reached the end of the chromosome
				if (nextRun == NULL) {
					break;
				}

				// Check for self-cyclic run
				if (nextRun == prevRun) {
					cerr << "Self-cyclic run found!" << endl;
					break;
				}

				// Check if the next run is the first run on the
				// chromosome, indicating a cyclic path.  If there is
				// a cyclic path we'll allow one more join attempt and
				// then break.
				bool cyclicPath = nextRun == chrom->getFirstRun();
				
				// Match the next run's left side with the prev run's
				// right side
				if (nextRun->getLeftChrom(g) != chrom or
					not nextRun->isLeftForward(g)) {
					nextRun->flip();
					assert(nextRun->getLeftChrom(g) == chrom and
						   nextRun->isLeftForward(g));
				}
				
// 				cerr << "Attempting to join runs:" << '\n'
// 					 << "Prev run:" << '\n'
// 					 << *prevRun << '\n'
// 					 << "Mask: " << prevRun->getMask() << '\n'
// 					 << "Next run:" << '\n'
// 					 << *nextRun << '\n'
// 					 << "Mask: " << nextRun->getMask() << '\n';
				
				// Attempt to join next run to previous run if they
				// are adjacent in at least minAdjacent genomes
				if (prevRun->isLeftOf(*nextRun, minAdjacent)) {
					Run* joinedRun = nextRun->joinTo(*prevRun, maxDist);
					if (joinedRun != NULL) {
						joinedRun->claimCliques();
						nextRun = joinedRun;
// 						cerr << "Yes, these runs are joined" << '\n';
					} else {
// 						cerr << "No, we can not join these runs" << '\n';
					}
				}
				
				prevRun = nextRun;

				// Make sure we are still on the same chromosome
				if (prevRun->getRightChrom(g) != chrom) {
					break;
				}

				// Stop if we are travelling back along a cyclic path
				if (cyclicPath) {
					break;
				}
			}
		}
	}
}

void filterIntraRunEdges() {
	vector<Run*> runs;
	getRuns(runs);

	vector<Run*>::iterator it;
	for (it = runs.begin(); it != runs.end(); ++it) {
		if ((*it)->isSignificant()) {
			(*it)->filterIntraRunEdges();
		}
	}
}

void filterInterRunEdges() {
	vector<Run*> runs;
	getRuns(runs);
	for_each(runs.begin(), runs.end(), mem_fun(&Run::filterInterRunEdges));
	for_each(runs.begin(), runs.end(),
			 bind2nd(mem_fun(&Run::setVisited), false));
}

void breakRuns() {
	vector<Run*> runs;
	getRuns(runs);

	vector<Run*>::iterator it;
	for (it = runs.begin(); it != runs.end(); ++it) {
		delete (*it);
	}
}

void removeInsignificantRuns() {
	vector<Run*> runs;
	getRuns(runs);

	vector<Run*>::iterator it;
	for (it = runs.begin(); it != runs.end(); ++it) {
		if (!(*it)->isSignificant()) {
			(*it)->deleteCliques();
			delete (*it);
		}
	}
}

void removeSingletons() {
	vector<Run*> runs;
	getRuns(runs);

	vector<Run*>::iterator it;
	for (it = runs.begin(); it != runs.end(); ++it) {
		(*it)->removeSingletons();
		if ((*it)->getLength() == 0) {
			delete (*it);
		}
	}

}	


void breakCyclicRuns() {
	// Check for cycles in each genome (only draft genomes, because
	// genomes not assembled by this program can not form cycles)
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		Genome* gen = Genome::getGenome(g);

		// Only consider draft genomes
		if (!gen->isDraft()) {
			continue;
		}

		// Make a mask with just this genome's flag set
		Mask included = 0;
		included[g] = 1;
				 
		// Get all of the runs that contain this genome
		vector<Run*> runs;
		getRuns(runs, included);

		// Examine every run
		vector<Run*>::iterator it;
		for (it = runs.begin(); it != runs.end(); ++it) {
			Run* start = *it;

			// Check if we have already examined this run
			if (start->wasVisited()) {
				continue;
			}

			start->setVisited();
			bool cycleFound = false;
			Run* lastAssembled = NULL; // The last run found during
									   // the traversal that was
									   // assembled.  This run will be
									   // split if a cycle is found.
			
			// Traverse twice, once to the right, and once to the left
			for (size_t i = 0; i < 2; ++i) {
				// Reverse the start run so that we go the opposite
				// direction on the second traversal
				if (i == 1) {
					start->flip();
				}
				
				Run* curr = start;
				while (true) {
					// Check if this run is assembled
					if (curr->isAssembled(g)) {
						lastAssembled = curr;
					}
					
					// Traverse to the right
					Run* next = curr->rightRun(g);

					// If we can not go any further to the right, stop
					if (next == NULL) {
						break;	
					}
					// If the next run has already been visited, it
					// must be the start run, and we must also be on
					// the first traversal
					else if (next->wasVisited()) {
						assert(i == 0 && next == start);
						cycleFound = true;
						break;
					}
					
					// Orient next run so that its left side matches
					// up with the right side of the current run
					if ((curr->getRightChrom(g) != next->getLeftChrom(g)) ||
						(curr->isRightForward(g) != next->isLeftForward(g))) {
						next->flip();
					}
					
					curr = next;
					curr->setVisited();
				}

				if (cycleFound) {
					break;
				}
			}

			// If a cycle was found, split off the right end of the
			// last assembled run
			if (cycleFound) {
				cerr << "Found cyclic run!" << '\n';
				assert(lastAssembled != NULL);
				Run* splitOff = lastAssembled->splitRunRight(g);
				splitOff->setVisited();
			}
		}

		// Re-find runs in this genome
		runs.clear();
		getRuns(runs, included);

		// Clear visited flags
		for (it = runs.begin(); it != runs.end(); ++it) {
			assert((*it)->wasVisited());
			(*it)->setVisited(false);
		}
	}
}

void Run::removeSingletons() {
	// Loop until we are no longer removing genomes or cliques
	bool cliquesRemoved = true;
	while (cliquesRemoved) {
		cliquesRemoved = false;

		// Remove singleton genomes
		bool genomeRemoved = false;
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if (hasGenome(g) && starts[g] == ends[g]) {
				starts[g]->removeAnchor(g);
				genomeRemoved = true;
			}
		}

		// If no genomes were removed, then no cliques will need to be removed
		if (!genomeRemoved) {
			break;
		}

		// Make a copy of the cliques originally in this run and reset
		vector<Clique*> oldCliques(cliques.begin(), cliques.end());
		reset();

		// Delete cliques of size 1 and add back in cliques of size > 1
		for (vector<Clique*>::iterator it = oldCliques.begin();
			 it != oldCliques.end(); ++it) {
			if ((*it)->getSize() >= 2 && (*it)->isConnected()) {
				addCliqueToRight(*it);
			} else {
				delete *it;
				cliquesRemoved = true;
			}
		}
	}
}

// Number the runs in RUNS in order starting at 1
void numberRuns(vector<Run*>& runs) {
	for (size_t i = 0; i < runs.size(); ++i) {
		runs[i]->setNum(i + 1);
	}
}
