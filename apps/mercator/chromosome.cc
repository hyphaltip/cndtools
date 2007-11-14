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

#include "chromosome.hh"
#include "anchor.hh"
#include "clique.hh"
#include "genome.hh"
#include "run.hh"

#include "bio/formats/agp.hh"

Chromosome::Chromosome(Genome* genome,
					   const string& name,
					   const GenomicDist length) :
	genome(genome),
	name(name),
	length(length),
	reversed(false)
{}

Genome* Chromosome::getGenome() const {
	return genome;
}

const string& Chromosome::getName() const {
	return name;
}

string Chromosome::getUniqueName() const {
	return genome->getName() + "_" + name;
}

GenomicDist Chromosome::getLength() const {
	return length;
}

GenomicDist Chromosome::getGenomeStart() const {
	return genomeStart;
}

GenomicDist Chromosome::getCoverage() const {
	GenomicDist covered = 0;
	for (Run* curr = getFirstRun(); curr != NULL;
		 curr = curr->nextRun(getGenome()->getNum())) {
		covered += curr->getCoverage(getGenome()->getNum());
	}
	return covered;
}

float Chromosome::getPctCoverage() const {
	return static_cast<float>(getCoverage()) / length;
}

bool Chromosome::isReversed() const {
	return reversed;
}

bool Chromosome::isPartOfAssembled() const {
	return !anchors.empty() && getFirstAnchor()->getChrom() != this;
}

void Chromosome::setGenomeStart(GenomicDist start) {
	genomeStart = start;
}

void Chromosome::addAnchor(Anchor* a) {
	anchors.push_back(a);
}

void Chromosome::initAnchors() {
	sort(anchors.begin(), anchors.end(), AnchorSorter());
	
	// Initialize next and prev pointers of anchors
	Anchor* prev = NULL;
	for (size_t i = 0; i < anchors.size(); ++i) {
		Anchor* next = anchors[i];
		next->setPrev(prev);
		next->setNext(NULL);
		if (prev != NULL) {
			prev->setNext(next);
		}
		prev = next;
	}
}

size_t Chromosome::getNumAnchors() const {
	return anchors.size();
}

size_t Chromosome::getNumAnchorsInCliques() const {
	return getNumCliques();
}

size_t Chromosome::getNumAnchorsInRuns() const {
	size_t numAnchors = 0;
	for (Run* curr = getFirstRun();
		 curr != NULL;
		 curr= curr->nextRun(getGenome()->getNum())) {
		numAnchors += curr->getNumAnchors(getGenome()->getNum());
	}
	return numAnchors;
}

size_t Chromosome::getNumCliques() const {
	size_t numCliques = 0;
	for (Clique* curr = getFirstClique();
		 curr != NULL;
		 curr = curr->nextClique(genome->getNum())) {
		++numCliques;
	}
	return numCliques;
}

size_t Chromosome::getNumRuns() const {
	size_t numRuns = 0;
	for (Run* curr = getFirstRun();
		 curr != NULL;
		 curr = curr->nextRun(genome->getNum())) {
		++numRuns;
	}
	return numRuns;
}

Anchor* Chromosome::getFirstAnchor() const {
	return (anchors.empty() ? NULL : anchors.front());
}

Anchor* Chromosome::getLastAnchor() const {
	return (anchors.empty() ? NULL : anchors.back());
}

Clique* Chromosome::getFirstClique() const {
	return (anchors.empty() ? NULL :
			(anchors.front()->isInClique() ?
			 anchors.front()->getClique() :
			 anchors.front()->nextClique()));
}

Clique* Chromosome::getLastClique() const {
	return (anchors.empty() ? NULL :
			(anchors.back()->isInClique() ?
			 anchors.back()->getClique() :
			 anchors.back()->nextClique(false)));
}

Run* Chromosome::getFirstRun() const {
	return (anchors.empty() ? NULL :
			(anchors.front()->isInRun() ?
			 anchors.front()->getRun() :
			 anchors.front()->nextRun(true)));
}

Run* Chromosome::getLastRun() const {
	return (anchors.empty() ? NULL :
			(anchors.back()->isInRun() ?
			 anchors.back()->getRun() :
			 anchors.back()->nextRun(false)));
}

void Chromosome::reverse() {
	reversed = !reversed;
	for_each(anchors.begin(), anchors.end(), mem_fun(&Anchor::reverse));
	std::reverse(anchors.begin(), anchors.end());
}

void Chromosome::unflipAnchors() {
	for_each(anchors.begin(), anchors.end(), mem_fun(&Anchor::unflip));
}

size_t Chromosome::getNumAnchorsRepetitive() const {
	return std::count_if(anchors.begin(), anchors.end(),
						 std::mem_fun(&Anchor::isMarked));
}

void Chromosome::writeRunPerm(std::ostream& strm) const {
	for (Run* curr = getFirstRun(); curr != NULL;
		 curr = curr->nextRun(genome->getNum())) {
		if (curr->isForward(genome->getNum())) {
			strm << curr->getNum() << ' ';
		} else {
			strm << -1 * curr->getNum() << ' ';
		}
	}
	strm << "$ # " << name << '\n';
}

void Chromosome::writeAGPLines(std::ostream& stream, size_t& recNum) const {
	stream << bio::formats::agp::Record(name, "", 1, length, recNum,  'D', 
										name, 1, length, '+');
}
