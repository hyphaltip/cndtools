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

#include "bio/phylogenetic/Tree.hh"
#include "bio/phylogenetic/TreeVisitor.hh"
#include "bio/phylogenetic/ConstTreeVisitor.hh"

#include <iostream>
#include <queue>
#include <limits>

namespace bio { namespace phylogenetic {

	Tree::Tree(const std::string& label, double branchLength)
		: label(label),
		  num(0),
		  branchLength(branchLength),
		  descendants() {
	}

	Tree::~Tree() {
		for (DescendantList::iterator it = descendants.begin();
			 it != descendants.end(); ++it) {
			delete *it;
		}
	}
	
	void Tree::dfsTraverse(ConstTreeVisitor& visitor) const {
		visitor.preVisitTree(this);
		for (DescendantList::const_iterator it = descendants.begin();
			 it != descendants.end(); ++it) {
			visitor.preVisitEdge(this, *it);
			(*it)->dfsTraverse(visitor);
			visitor.postVisitEdge(this, *it);
		}
		visitor.postVisitTree(this);
	}

	void Tree::dfsTraverse(TreeVisitor& visitor) {
		visitor.preVisitTree(this);
		for (DescendantList::iterator it = descendants.begin();
			 it != descendants.end(); ++it) {
			visitor.preVisitEdge(this, *it);
			(*it)->dfsTraverse(visitor);
			visitor.postVisitEdge(this, *it);
		}
		visitor.postVisitTree(this);
	}
	
	void Tree::bfsTraverse(ConstTreeVisitor& visitor) const {
		std::queue<const Tree*> unprocessed;
		unprocessed.push(this);
		while (!unprocessed.empty()) {
			const Tree* next = unprocessed.front();
			visitor.preVisitTree(next);
			for (DescendantList::const_iterator it = next->descendants.begin();
				 it != next->descendants.end(); ++it) {
				unprocessed.push(*it);
				visitor.postVisitEdge(next, *it);				
			}
			unprocessed.pop();
			visitor.postVisitTree(next);
		}
	}

	void Tree::bfsTraverse(TreeVisitor& visitor) {
		std::queue<Tree*> unprocessed;
		unprocessed.push(this);
		while (!unprocessed.empty()) {
			Tree* next = unprocessed.front();
			visitor.preVisitTree(next);
			for (DescendantList::iterator it = next->descendants.begin();
				 it != next->descendants.end(); ++it) {
				unprocessed.push(*it);
				visitor.postVisitEdge(next, *it);				
			}
			unprocessed.pop();
			visitor.postVisitTree(next);
		}
	}

	Tree* Tree::copy() const {
		Tree* t = new Tree(label, branchLength);
		for (size_t i = 0; i < descendants.size(); ++i) {
			t->addDescendant(descendants[i]->copy());
		}
		return t;
	}
	
	Tree* Tree::getSubtree(util::stl::hash_set<std::string>& taxa,
						   bool keepDescendants) const {
		bool includedTaxon = (taxa.find(label) != taxa.end());
		
		if (keepDescendants and includedTaxon) {
			return copy();
		}

		DescendantList subtreeDescendants;
		for (size_t i = 0; i < descendants.size(); ++i) {
			Tree* t = descendants[i]->getSubtree(taxa, keepDescendants);
			if (t != NULL) {
				subtreeDescendants.push_back(t);
			}
		}

		if (includedTaxon or subtreeDescendants.size() > 1) {
			Tree* t = new Tree(label, branchLength);
			for (size_t i = 0; i < subtreeDescendants.size(); ++i) {
				t->addDescendant(subtreeDescendants[i]);
			}
			return t;
		} else if (subtreeDescendants.empty()) {
			return NULL;
		} else {
			Tree* child = subtreeDescendants.front();
			child->setBranchLength(child->getBranchLength() + branchLength);
			return child;
		}
	}
	
	class TreeNumbererVisitor : public TreeVisitor {
	private:
		unsigned int counter;
	public:
		TreeNumbererVisitor() : counter(0) {}
		void postVisitTree(Tree* t) { t->setNum(counter++); }
	};

	void Tree::setNumbers() {
		TreeNumbererVisitor numberer;
		dfsTraverse(numberer);
	}

    size_t Tree::getNumLeaves() const {
		if (isLeaf()) {
			return 1;
		} else {
			size_t numLeaves = 0;
			DescendantList::const_iterator it;
			for (it = descendants.begin(); it != descendants.end(); ++it) {
				numLeaves += (*it)->getNumLeaves();
			}
			return numLeaves;
		}
	}
		
    size_t Tree::getNumNodes() const {
		size_t numNodes = 1;
		DescendantList::const_iterator it;	
		for (it = descendants.begin(); it != descendants.end(); ++it) {
			numNodes += (*it)->getNumNodes();
		}
		return numNodes;
	}
	

	class TreePairwiseDistanceVisitor : public ConstTreeVisitor {
	private:
		util::Matrix<double>& m;
	public:
		TreePairwiseDistanceVisitor(util::Matrix<double>& m)
			: m(m)
		{
			for (size_t i = 0; i < m.getNumRows(); ++i) {
				for (size_t j = 0; j < m.getNumCols(); ++j) {
					m(i, j) = (i == j ? 0.0 : std::numeric_limits<double>::min());
				}
			}
		}

		void preVisitEdge(const Tree* parent, const Tree* child) {
			for (size_t i = 0; i < m.getNumCols(); ++i) {
				if (i == child->getNum()) {
					continue;
				}
				double parentDistance = m(parent->getNum(), i);
				if (parentDistance != std::numeric_limits<double>::min()) {
					double childDistance = parentDistance + child->getBranchLength();
					m(child->getNum(), i) = childDistance;
					m(i, child->getNum()) = childDistance;					
				}
			}
		}
	};
	
	void Tree::getPairwiseDistances(util::Matrix<double>& distances) const {
		size_t numNodes = getNumNodes();
		distances.resize(numNodes, numNodes);
		TreePairwiseDistanceVisitor visitor(distances);
		dfsTraverse(visitor);
	}

	class TreeOrdererVisitor : public ConstTreeVisitor {
	private:
		std::vector<const Tree*>& order;
	public:
		TreeOrdererVisitor(std::vector<const Tree*>& order) :  order(order) {}
		void postVisitTree(const Tree* t) { order.push_back(t); }
	};
	
	void Tree::getTopologicalOrder(std::vector<const Tree*>& order) const {
		TreeOrdererVisitor visitor(order);
		dfsTraverse(visitor);
	}

} }
