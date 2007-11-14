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

#ifndef __EDGE_HH__
#define __EDGE_HH__

#include "types.hh"
#include "anchor.hh"

class Edge {
public:
	typedef int ScoreType;
	
	Edge(Anchor* a1,
		 Anchor* a2,
		 const ScoreType score);
	
	Anchor* getAnchor1() const;
	Anchor* getAnchor2() const;
	Anchor* getOtherAnchor(const Anchor* a) const;
	size_t getOtherGenomeNum(const Anchor* a) const;
	bool hasGenomes(const size_t g1, const size_t g2) const;
	bool isActive() const;
	ScoreType getScore() const;

	bool shouldBePruned(const float prunePct) const;
	
	void addEdge();

	friend ostream& operator<<(ostream& strm, const Edge& e);
	
	void print();

private:
	Anchor* a1;
	Anchor* a2;
	ScoreType score;
};

struct EdgeSorter {
	bool operator()(const Edge* e1, const Edge* e2) const {
		return e1->getScore() < e2->getScore();
	}
};

inline Anchor* Edge::getAnchor1() const { return a1; }

inline Anchor* Edge::getAnchor2() const { return a2; }

inline Anchor* Edge::getOtherAnchor(const Anchor* a) const {
	return a1 == a ? a2 : a1;
}

inline size_t Edge::getOtherGenomeNum(const Anchor* a) const {
	return getOtherAnchor(a)->getGenomeNum();
}

inline bool Edge::isActive() const {
	return a1->hasEdge(this);
}

inline Edge::ScoreType Edge::getScore() const { return score; }

inline void Edge::addEdge() {
	a1->addEdge(this);
	a2->addEdge(this);
}

#endif // __EDGE_HH__
