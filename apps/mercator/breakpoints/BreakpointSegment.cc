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

#include <iostream>
#include <cassert>

#include "BreakpointSegment.hh"
#include "Segment.hh"
#include "Genome.hh"

size_t BreakpointSegment::resolution;

void BreakpointSegment::setResolution(size_t resolution) {
	BreakpointSegment::resolution = resolution;
}

size_t BreakpointSegment::getResolution() {
	return resolution;
}

BreakpointSegment::BreakpointSegment(Segment* upstreamSegment,
									 Segment* downstreamSegment)
	: upstreamSegment(upstreamSegment),
	  downstreamSegment(downstreamSegment),
	  changed(true)
{
	if (upstreamSegment->getStrand() == '+') {
		upstreamSegment->setDownstreamSegment(this);
	} else {
		upstreamSegment->setUpstreamSegment(this);
	}
	if (downstreamSegment->getStrand() == '+') {
		downstreamSegment->setUpstreamSegment(this);
	} else {
		downstreamSegment->setDownstreamSegment(this);
	}
}

bool BreakpointSegment::isColinearWith(const BreakpointSegment& other) const {
	return ((upstreamSegment->getSegSet() == other.upstreamSegment->getSegSet() and
			 downstreamSegment->getSegSet() == other.downstreamSegment->getSegSet()) or
			(upstreamSegment->getSegSet() == other.downstreamSegment->getSegSet() and
			 downstreamSegment->getSegSet() == other.upstreamSegment->getSegSet()));
}

Genome* BreakpointSegment::getGenome() const {
	return upstreamSegment->getGenome();
}

std::string BreakpointSegment::getChrom() const {
	return upstreamSegment->getChrom();
}

size_t BreakpointSegment::getStart() const {
	return upstreamSegment->getEnd();
}	

size_t BreakpointSegment::getEnd() const {
	return downstreamSegment->getStart();
}

bool BreakpointSegment::isChanged() const {
	return changed;
}

std::string BreakpointSegment::getSeq(char strand) const {
	return getGenome()->getSeq(getChrom(), getStart(), getEnd(), strand);
}

void BreakpointSegment::setNum(size_t num) { this->num = num; }
size_t BreakpointSegment::getNum() const { return num; }

size_t BreakpointSegment::getLength() const { return getEnd()  - getStart(); }

size_t BreakpointSegment::getBoundedLength() const {
	return upperBound - lowerBound;
}

void BreakpointSegment::resetScores() {
	scores.resize(std::min(upperBound - lowerBound + 1, resolution));
	std::fill(scores.begin(), scores.end(), 0);
}

void BreakpointSegment::initBounds() {
	lowerBound = 0;
	upperBound = getLength();
}

void BreakpointSegment::updateBounds(size_t minBoundedLength) {
	if (getBoundedLength() <= minBoundedLength) {
		return;
	}

	size_t lowHalfLength = minBoundedLength / 2;
	size_t upHalfLength = minBoundedLength / 2 + (minBoundedLength % 2);
	size_t midpoint = getPosition(bestIndex);
	
	if ((midpoint - lowerBound) < lowHalfLength) {
		upperBound = lowerBound + minBoundedLength;
	} else if ((upperBound - midpoint) < upHalfLength) {
		lowerBound = upperBound - minBoundedLength;
	} else {
		lowerBound = midpoint - lowHalfLength;
		upperBound = midpoint + upHalfLength;
	}

	assert(upperBound <= getLength());
}

// void BreakpointSegment::updateBounds() {
// 	if (scores.size() < resolution) {
// 		return;
// 	}

// 	size_t newLowerBound = lowerBound;
// 	size_t newUpperBound = upperBound;
// 	if (bestIndex == 0) {
// 		newUpperBound = getPosition(2);
// 	} else if (bestIndex == scores.size() - 1) {
// 		newLowerBound = getPosition(scores.size() - 3);
// 	} else {
// 		newLowerBound = getPosition(bestIndex - 1);
// 		newUpperBound = getPosition(bestIndex + 1);
// 	}

// 	lowerBound = newLowerBound;
// 	upperBound = newUpperBound;
// }

void BreakpointSegment::setBestBreakpoint() {
	assert(scores.size() < resolution);
	lowerBound = lowerBound + bestIndex;
	upperBound = lowerBound;
}

void BreakpointSegment::setRelBreakpoint(size_t bp) {
	lowerBound = bp;
	upperBound = lowerBound;
}

size_t BreakpointSegment::getRelBreakpoint() const {
	return lowerBound;
}

size_t BreakpointSegment::getBreakpoint() const {
	return getStart() + lowerBound;
}

size_t BreakpointSegment::getNumIndices() const {
	return scores.size();
}

size_t BreakpointSegment::getPosition(unsigned long long i) const {
	if (i == 0) {
		return lowerBound;
	} else {
		return lowerBound + (i * (upperBound - lowerBound)) / (scores.size() - 1);
	}
}

Score BreakpointSegment::getScore(size_t i) const {
	try {
		return scores.at(i);
	} catch (std::runtime_error& e) {
		throw std::runtime_error(std::string("getScore: ") + e.what());
	}
}

void BreakpointSegment::setScore(size_t i, Score score) {
	try {
		scores.at(i) = score;
	} catch (std::runtime_error& e) {
		throw std::runtime_error(std::string("getScore: ") + e.what());
	}
}

size_t BreakpointSegment::getBestIndex() const {
	return bestIndex;
}

void BreakpointSegment::setBestIndex(size_t i) {
	bestIndex = i;
}

void BreakpointSegment::split() {
	//std::cerr << "Splitting: " << *this;
	size_t breakpoint = getBreakpoint();
	upstreamSegment->setEnd(breakpoint);
	downstreamSegment->setStart(breakpoint);

	//std::cerr << " breakpoint = " << breakpoint << '\n';
}

std::ostream& operator<<(std::ostream& strm, const BreakpointSegment& seg) {
	return strm << '('
				<< seg.getGenome()->getName() << ':'
				<< seg.getChrom() << ':'
				<< seg.getStart() << '-'
				<< seg.getEnd()
				<< ')';
}
