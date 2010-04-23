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

#ifndef __BIO_PHYLOGENETIC_TREE_HH__
#define __BIO_PHYLOGENETIC_TREE_HH__

#include <string>
#include <vector>

#include "boost/unordered_set.hpp"
#include "util/matrix.hh"
		
namespace bio { namespace phylogenetic {

	class TreeVisitor;
	class ConstTreeVisitor;

	// A node in a phylogenetic tree
	class Tree {
	private:
		typedef std::vector<Tree*> DescendantList;
		std::string label;
		size_t num;
		double branchLength;
		DescendantList descendants;

	public:
		// Construct from optional node label and optional branch
		// length from parent node
		Tree(const std::string& label = "", double branchLength = 0.0);

		// Destruct
		~Tree();

		// Returns true if this node is a leaf node
		bool isLeaf() const { return descendants.empty(); }

		// Add tree D as a descendant of this tree node
		void addDescendant(Tree* d) { descendants.push_back(d); }

		// Get the label of this node
		const std::string& getLabel() const { return label; }

		// Get the length of the branch connecting this node to its parent
		double getBranchLength() const { return branchLength; }

		// Get the topological number of this node.  The method
		// setNumbers must have been called previously for this number
		// to be valid.
		size_t getNum() const { return num; }

		// Get the number of nodes in this (sub)tree
		size_t getNumNodes() const;

		// Get the number of leaf nodes in this (sub)tree
		size_t getNumLeaves() const;

		// Get the number of descendants of this tree node
		size_t getNumDescendants() const { return descendants.size(); }

		// Get the ith descendant
		Tree* getDescendant(size_t i) const { return descendants.at(i); }
		
		// Get the descendants of this tree node
		const DescendantList& getDescendants() const { return descendants; }

		// Set the label of this node
		void setLabel(const std::string& l) { label = l; }

		// Set the length of the branch connecting this node to its parent
		void setBranchLength(double bl) { branchLength = bl; }

		// Set the topological number of this node.  Used by
		// setNumbers method
		void setNum(size_t n) { num = n; }

		// Traverse this tree in depth-first order with VISITOR
		void dfsTraverse(TreeVisitor& visitor);

		// Traverse this tree in breadth-first order with VISITOR
		void bfsTraverse(TreeVisitor& visitor);		

		// Traverse this tree in depth-first order with VISITOR.
		// VISITOR may not modify the tree.
		void dfsTraverse(ConstTreeVisitor& visitor) const;

		// Traverse this tree in breadth-first order with VISITOR.
		// VISITOR may not modify the tree.
		void bfsTraverse(ConstTreeVisitor& visitor) const;

		// Set the topological numbers of the nodes in this tree
		void setNumbers();

		// Get the pairwise distances between all nodes in the matrix
		// DISTANCES
		void getPairwiseDistances(util::Matrix<double>& distances) const;

		// Get all of the nodes in this tree in toplogical order in
		// the vector ORDER
		void getTopologicalOrder(std::vector<const Tree*>& order) const;

		Tree* copy() const;
		
		// Return tree with given TAXA contained within this tree
		Tree* getSubtree(boost::unordered_set<std::string>& taxa,
						 bool keepDescendants = false) const;
	};

} }

#endif // __BIO_PHYLOGENETIC_TREE_HH__
