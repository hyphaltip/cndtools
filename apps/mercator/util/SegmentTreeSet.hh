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

#include <vector>
#include <map>

#include "bio/genome/BasicInterval.hh"
#include "boost/unordered_map.hpp"

template<typename Map1, typename Map2>
void invert_map(const Map1& m1, Map2& m2) {
	typedef typename Map1::const_iterator M1Iter;
	for (M1Iter m1_iter = m1.begin(); m1_iter != m1.end(); ++m1_iter) {
		m2[m1_iter->second] = m1_iter->first;
	}
}

class SegmentTreeSet {
public:
	class Sequence;
	class TreeNode;
	class Ancestor;
	class Leaf;

	typedef size_t Position;
	typedef boost::unordered_map<std::string, Sequence*> SeqMap;
	typedef SeqMap::iterator SeqMapIter;
	typedef SeqMap::const_iterator SeqMapConstIter;
	typedef std::vector<Sequence*> SeqList;
	typedef SeqList::iterator SeqListIter;
	typedef SeqList::const_iterator SeqListConstIter;

	typedef std::pair<TreeNode*, bool> TreeOrientPair;
	typedef std::pair<TreeNode*, TreeNode*> TreeNodePair;
	typedef std::map<Position, Leaf*> BlockMap;
	typedef BlockMap::value_type BlockMapValue;
	typedef BlockMap::iterator BlockMapIter;
	typedef BlockMap::const_iterator BlockMapConstIter;
	typedef std::pair<BlockMapIter, BlockMapIter> BlockMapIterPair;
	typedef boost::unordered_set<TreeNode*> TreeSet;
	typedef std::vector<Leaf*> SegmentList;

	SegmentTreeSet() {}

	// Add a sequence of name NAME and length LENGTH to the set
	Sequence* addSeq(const std::string& name, Position length);

	// Returns the number of sequences in the set
	size_t getNumSeqs() const { return seqList.size(); }

	// Returns the Ith sequence of the set
	Sequence* getSeq(size_t i) { return seqList[i]; }

	// Returns the sequence named NAME in the set
	Sequence* getSeq(const std::string& name) {
		assert(hasSeq(name));
		return seqMap[name];
	}

	// Returns true if the set has a sequence of name NAME
	bool hasSeq(const std::string& name) {
		return seqMap.find(name) != seqMap.end();
	}

	// Label the trees in the set from 1 to the number of trees
    void getTrees(TreeSet& trees) const;

	// Label the segments in the set
	void getSegments(SegmentList& segments) const;

	// Add a match of length LENGTH between S1 and S2, starting at
	// positions START1 and START2, respectively.  Coordinates are
	// relative to the reverse strand if FORWARD is set to false.
	void addMatch(Sequence& s1, Position start1, bool forward1,
				  Sequence& s2, Position start2, bool forward2,
				  Position length);

	// Add a match of length LENGTH between S1 and S2, starting at
	// positions START1 and START2, respectively.  Coordinates are
	// relative to the reverse strand if FORWARD is set to false.
	void addMatch(const std::string& s1, Position start1, bool forward1,
				  const std::string& s2, Position start2, bool forward2,
				  Position length);
	
	class Sequence {
	public:
		Sequence(const std::string& name, Position length);

		BlockMapIter cut(Position p);
		BlockMapIterPair cut(Position p1, Position p2, bool forward);

		bool isRangeRepetitive(Position p1,
							   Position p2,
							   bool forward,
							   size_t max_leaves);
		
		void getTrees(TreeSet& trees) const;
		void getSegments(SegmentList& segments) const;

		std::string name;
		Position length;
		BlockMap blocks;

	private:
		Position& flip(Position& p) const;
		void flip(Position& p1, Position& p2) const;
		void checkPosition(Position p) const;
	};

	class TreeNode {
	public:
		TreeNode() : parent(NULL), forward(true) {}
		TreeNode(bool forward) : parent(NULL), forward(forward) {}
		virtual ~TreeNode() {}
		virtual TreeNode* getChild1() const = 0;
		virtual TreeNode* getChild2() const = 0;
		virtual bool isLeaf() const = 0;

		TreeNode* getParent() const { return parent; }
		void setParent(TreeNode* parent) { this->parent = parent; }
		bool isRoot() const { return parent == NULL; }

		TreeNode* getRoot() const;
		std::pair<TreeNode*, bool> getRootAndForward() const;

		void flip() { forward = not forward; }

		virtual size_t getNumLeaves() const = 0;
		
		virtual size_t getLength() const = 0;
		virtual void orientTree() = 0;
		virtual std::string toString() const = 0;
		virtual TreeNodePair splitTree(const Position offset) = 0;
		virtual void getLeaves(std::vector<Leaf*>& leaves) const = 0;

		virtual bool hasNeighboringPositions() = 0;
		
	protected:
		TreeNode* parent;
		bool forward;
	};

	class Leaf : public TreeNode {
	public:
		Leaf(Sequence& seq, BlockMapIter iter, Position pos, bool forward);
		TreeNode* getChild1() const { return NULL; }
		TreeNode* getChild2() const { return NULL; }
		bool isLeaf() const { return true; }
		void split(const Position pos);
		TreeNode* join(Leaf* other, bool same_orientation);
		bool treeHasDuplicates() const;

		void orientTree() { return; }
		TreeNodePair splitTree(const Position offset);
		std::string toString() const;

		size_t getLength() const;

		size_t getNumLeaves() const { return 1; }
		
		void getLeaves(std::vector<Leaf*>& leaves) const;
		Position getStart() const { return iter->first; }
		bio::genome::BasicInterval getInterval() const;
		bool isUnique() const { return parent == NULL; }
		bool isForward() const { return forward; }

		const std::string& getSeqName() const { return seq.name; }

		bool inSameTreeAsNeighbor();
		bool hasNeighboringPositions() { return inSameTreeAsNeighbor(); }
		
		Leaf* nextLeaf();
		TreeNode* nextTree();
		Leaf* prevLeaf();
		TreeNode* prevTree();
	
	private:
		Position getPosition(const Position offset) const;
		Position getOffset(const Position pos) const;
	
		Sequence& seq;
		BlockMapIter iter;
	};

	class Ancestor : public TreeNode {
	public:
		Ancestor(TreeNode* child1, TreeNode* child2);
		TreeNode* getChild1() const { return child1; }
		TreeNode* getChild2() const { return child2; }
		bool isLeaf() const { return false; }

		void orientTree();
		TreeNodePair splitTree(const Position offset);
		std::string toString() const;
		size_t getNumLeaves() const;
		void getLeaves(std::vector<Leaf*>& leaves) const;
		size_t getLength() const;

		bool hasNeighboringPositions() {
			return child1->hasNeighboringPositions() or
				child2->hasNeighboringPositions();
		}
		
	private:
		TreeNode* child1;
		TreeNode* child2;
	};
	
protected:
	SeqList seqList;
	SeqMap seqMap;

	static void joinBlocks(BlockMapIterPair& bounds1,
						   BlockMapIterPair& bounds2,
						   bool forward1,
						   bool forward2);

	static BlockMapIter next_iter(BlockMapIter i, bool forward);
	static Leaf* next_block(BlockMapIter i, bool forward);
	static Position distance(BlockMapIter i1,
							 BlockMapIter i2,
							 bool forward);
	static Position position(BlockMapIter i, Position offset, bool forward);
	static void print_positions(BlockMapIterPair& bounds, bool forward);
	
};
