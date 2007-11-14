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

#include "edge.hh"
#include "anchor.hh"

Edge::Edge(Anchor* a1,
		   Anchor* a2,
		   const ScoreType score)
	: a1(a1),
	  a2(a2),
	  score(score) {
}

bool Edge::shouldBePruned(const float prunePct) const {
	return (a1->hasEdgesTo(a2->getGenomeNum()) &&
			score < prunePct * a1->getBestEdge(a2->getGenomeNum())->getScore())
		|| (a2->hasEdgesTo(a1->getGenomeNum()) &&
			score < prunePct * a2->getBestEdge(a1->getGenomeNum())->getScore());
}

ostream& operator<<(ostream& strm, const Edge& e) {
	return strm << *(e.a1) << '\t' << *(e.a2) << '\t' << e.score << '\n';
}

void Edge::print() {
	std::cerr << (*a1) << '\t' << (*a2) << '\t' << score << '\n';
}

bool Edge::hasGenomes(const size_t g1,
					  const size_t g2) const {
	return (a1->getGenomeNum() == g1 && a2->getGenomeNum() == g2)
		|| (a1->getGenomeNum() == g2 && a2->getGenomeNum() == g1);
}
