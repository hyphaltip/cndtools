#include <iostream>

#include "util/string.hh"

#include "SegmentTreeSet.hh"

std::pair<SegmentTreeSet::TreeNode*, bool>
SegmentTreeSet::TreeNode::getRootAndForward() const {
	TreeNode* root = const_cast<TreeNode*>(this);
	bool orient = forward;
	while (root->parent != NULL) {
		root = root->parent;
		if (not root->forward) { orient = not orient; }
	}
	return std::make_pair(root, orient);
}

SegmentTreeSet::TreeNode*
SegmentTreeSet::TreeNode::getRoot() const {
	TreeNode* root = const_cast<TreeNode*>(this);
	while (root->parent != NULL) { root = root->parent; }
	return root;
}

bool SegmentTreeSet::Leaf::treeHasDuplicates() const {
	std::vector<Leaf*> leaves;
	getRoot()->getLeaves(leaves);
	std::sort(leaves.begin(), leaves.end());
	for (size_t i = 1; i < leaves.size(); ++i) {
		if (leaves[i] == leaves[i - 1]) {
			return true;
		}
	}
	return false;
}

size_t SegmentTreeSet::Leaf::getLength() const {
	return (++BlockMapConstIter(iter))->first - iter->first;
}

void SegmentTreeSet::Leaf::split(const Position pos) {
	TreeNode* root = getRoot();
	root->orientTree();
	root->splitTree(getOffset(pos));
}

SegmentTreeSet::Leaf::Leaf(Sequence& seq,
						   BlockMapIter iter,
						   Position pos,
						   bool forward)
	: TreeNode(forward),
	  seq(seq),
	  iter(seq.blocks.insert(iter, BlockMapValue(pos, this)))
{
	assert(this->iter->first == pos);
}

bool
SegmentTreeSet::Leaf::inSameTreeAsNeighbor() {
	TreeNode* tree = getRoot();
	return (tree == nextTree() or tree == prevTree());
}

SegmentTreeSet::Leaf*
SegmentTreeSet::Leaf::nextLeaf() {
	return (++BlockMapConstIter(iter))->second;
}

SegmentTreeSet::TreeNode*
SegmentTreeSet::Leaf::nextTree() {
	Leaf* next = nextLeaf();
	return (next ? next->getRoot() : NULL);
}

SegmentTreeSet::Leaf*
SegmentTreeSet::Leaf::prevLeaf() {
	return (iter->first == 0 ? NULL : (--BlockMapConstIter(iter))->second);
}

SegmentTreeSet::TreeNode*
SegmentTreeSet::Leaf::prevTree() {
	Leaf* prev = prevLeaf();
	return (prev ? prev->getRoot() : NULL);
}

SegmentTreeSet::TreeNode*
SegmentTreeSet::Leaf::join(Leaf* other, bool same_orientation) {
	TreeOrientPair t1 = this->getRootAndForward();
	TreeOrientPair t2 = other->getRootAndForward();
	
	bool curr_same_orientation = (t1.second == t2.second);
	if (t1.first == t2.first) {
		return (curr_same_orientation == same_orientation ? t1.first : NULL);
	} else {
		if (curr_same_orientation != same_orientation) {
			t2.first->flip();
		}
		return new Ancestor(t1.first, t2.first);
	}
}

SegmentTreeSet::TreeNodePair
SegmentTreeSet::Leaf::splitTree(const Position offset) {
	Position p = getPosition(offset);
	Leaf* right = new Leaf(seq, iter, p, forward);
	return (forward ? TreeNodePair(this, right) : TreeNodePair(right, this));
}

std::string SegmentTreeSet::Leaf::toString() const {
	return util::string::toString(getInterval());
}	

SegmentTreeSet::Position
SegmentTreeSet::Leaf::getPosition(const Position offset) const {
	if (forward) {
		return iter->first + offset;
	} else {
		return (++BlockMapIter(iter))->first - offset;
	}
}

SegmentTreeSet::Position
SegmentTreeSet::Leaf::getOffset(const Position pos) const {
	if (forward) {
		return pos - iter->first;
	} else {
		return (++BlockMapIter(iter))->first - pos;
	}
}

void SegmentTreeSet::Leaf::getLeaves(std::vector<Leaf*>& leaves) const {
	leaves.push_back(const_cast<Leaf*>(this));
}

bio::genome::BasicInterval SegmentTreeSet::Leaf::getInterval() const {
	return bio::genome::BasicInterval(seq.name,
									  iter->first,
									  iter->first + getLength(),
									  (forward ? '+' : '-'));
}

SegmentTreeSet::Ancestor::Ancestor(TreeNode* child1, TreeNode* child2)
	: TreeNode(),
	  child1(child1),
	  child2(child2)
{
	child1->setParent(this);
	child2->setParent(this);
}

void SegmentTreeSet::Ancestor::orientTree() {
	if (not forward) {
		child1->flip();
		child2->flip();
		forward = true;
	}
	child1->orientTree();
	child2->orientTree();
}

std::string SegmentTreeSet::Ancestor::toString() const {
	return "(" + child1->toString() + "," + child2->toString() + ")";
}

void SegmentTreeSet::Ancestor::getLeaves(std::vector<Leaf*>& leaves) const {
	child1->getLeaves(leaves);
	child2->getLeaves(leaves);
}

size_t SegmentTreeSet::Ancestor::getNumLeaves() const {
	return child1->getNumLeaves() + child2->getNumLeaves();
}

size_t SegmentTreeSet::Ancestor::getLength() const {
	return child1->getLength();
}

SegmentTreeSet::TreeNodePair
SegmentTreeSet::Ancestor::splitTree(const Position offset) {
	TreeNodePair trees1 = child1->splitTree(offset);
	TreeNodePair trees2 = child2->splitTree(offset);
	child1 = trees1.first;
	child1->setParent(this);
	child2 = trees2.first;
	child2->setParent(this);
	return TreeNodePair(this, new Ancestor(trees1.second, trees2.second));
}

SegmentTreeSet::Sequence*
SegmentTreeSet::addSeq(const std::string& name, Position length) {
	Sequence* s = new Sequence(name, length);
	seqList.push_back(s);
	seqMap[name] = s;
	return s;
}

void SegmentTreeSet::getTrees(TreeSet& trees) const {
	trees.clear();
	for (SeqListConstIter s = seqList.begin(); s != seqList.end(); ++s) {
		(*s)->getTrees(trees);
	}
}

void SegmentTreeSet::getSegments(SegmentList& segments) const {
	segments.clear();
	for (SeqListConstIter s = seqList.begin(); s != seqList.end(); ++s) {
		(*s)->getSegments(segments);
	}
}

void SegmentTreeSet::Sequence::getTrees(TreeSet& trees) const {
	for (BlockMapConstIter b = blocks.begin(); b != blocks.end(); ++b) {
		if (b->second != NULL) {
			trees.insert(b->second->getRoot());
		}
	}
}

void SegmentTreeSet::Sequence::getSegments(SegmentList& segments) const {
	for (BlockMapConstIter b = blocks.begin(); b != blocks.end(); ++b) {
		if (b->second != NULL) {
			segments.push_back(b->second);
		}
	}
}

SegmentTreeSet::Sequence::Sequence(const std::string& name, Position length)
	: name(name), length(length), blocks() {
	new Leaf(*this, blocks.begin(), 0, true);
	blocks.insert(BlockMapValue(length, NULL));
}

void SegmentTreeSet::Sequence::checkPosition(Position p) const {
	using util::string::toString;
	if (p > length) {
		throw std::runtime_error("Position greater than sequence length: "
								 "seq: " + name +
								 " length: " + toString(length) +
								 " position: " + toString(p));
	}
}

SegmentTreeSet::BlockMapIter
SegmentTreeSet::Sequence::cut(Position p) {
	checkPosition(p);

	BlockMapIter i = blocks.lower_bound(p);
	if (i->first != p) {
		(--i)->second->split(p);
		++i;
	}
	if (i->first != p) {
		std::cerr << p << ' '
				  << i->first << ' '
				  << (--BlockMapIter(i))->first << ' '
				  << (++BlockMapIter(i))->first << '\n';
	}
	assert(i->first == p);
	return i;
}


bool
SegmentTreeSet::Sequence::isRangeRepetitive(Position p1,
											Position p2,
											bool forward,
											size_t max_leaves) {
	checkPosition(p1);
	checkPosition(p2);
	
	if (not forward) { flip(p1); flip(p2); }
	BlockMapIter i = blocks.upper_bound(p1);
	if (i == blocks.end()) { return false; }
	--i;
	while (p2 > i->first) {
		if (i->second->getRoot()->getNumLeaves() > max_leaves) {
			return true;
		}
		++i;
	}
	return false;
}

SegmentTreeSet::BlockMapIterPair
SegmentTreeSet::Sequence::cut(Position p1, Position p2, bool forward) {
	if (not forward) { flip(p1); flip(p2); }
	return BlockMapIterPair(cut(p1), cut(p2));
}

SegmentTreeSet::Position&
SegmentTreeSet::Sequence::flip(Position& p) const {
	return p = length - p;
}

void SegmentTreeSet::Sequence::flip(Position& p1, Position& p2) const {
	std::swap(flip(p1), flip(p2));
}

inline
SegmentTreeSet::BlockMapIter
SegmentTreeSet::next_iter(BlockMapIter i, bool forward) {
	return (forward ? ++i : --i);
}

inline
SegmentTreeSet::Leaf*
SegmentTreeSet::next_block(BlockMapIter i, bool forward) {
	return (forward ? i->second : (--i)->second);
}

inline
SegmentTreeSet::Position
SegmentTreeSet::distance(BlockMapIter i1,
						 BlockMapIter i2,
						 bool forward) {
	return (forward ? i2->first - i1->first : i1->first - i2->first);
}

inline
SegmentTreeSet::Position
SegmentTreeSet::position(BlockMapIter i, Position offset, bool forward) {
	return (forward ? i->first + offset : i->first - offset);
}

void SegmentTreeSet::print_positions(BlockMapIterPair& bounds, bool forward) {
	SegmentTreeSet::BlockMapIter curr(bounds.first);
	while (curr != bounds.second) {
		std::cerr << curr->first << ' ';
		curr = next_iter(curr, forward);
	}
	std::cerr << '\n';
}

void SegmentTreeSet::addMatch(const std::string& s1,
							  Position start1,
							  bool forward1,
							  const std::string& s2,
							  Position start2,
							  bool forward2,
							  Position length) {
// 	std::cerr << "Adding match: "
// 			  << s1 << ' '
// 			  << start1 << ' '
// 			  << forward1 << ' '
// 			  << s2 << ' '
// 			  << start2 << ' '
// 			  << forward2 << ' '
// 			  << length << '\n';
	
	addMatch(*getSeq(s1), start1, forward1,
			 *getSeq(s2), start2, forward2,
			 length);
}

void SegmentTreeSet::addMatch(Sequence& s1, Position start1, bool forward1,
							  Sequence& s2, Position start2, bool forward2,
							  Position size) {
	if (s1.isRangeRepetitive(start1, start1 + size, forward1, 100) or
		s2.isRangeRepetitive(start2, start2 + size, forward2, 100)) {
		return;
	}
	
	BlockMapIterPair bounds1 = s1.cut(start1, start1 + size, forward1);
	BlockMapIterPair bounds2 = s2.cut(start2, start2 + size, forward2);
	joinBlocks(bounds1, bounds2, forward1, forward2);
}

void SegmentTreeSet::joinBlocks(BlockMapIterPair& bounds1,
								BlockMapIterPair& bounds2,
								bool forward1,
								bool forward2) {
	BlockMapIter curr1(bounds1.first), curr2(bounds2.first);
	while (curr1 != bounds1.second) {
		BlockMapIter next1(next_iter(curr1, forward1));
		BlockMapIter next2(next_iter(curr2, forward2));
		Position dist1 = distance(curr1, next1, forward1);
		Position dist2 = distance(curr2, next2, forward2);
		if (dist1 < dist2) {
			next_block(curr2, forward2)->split(position(curr2, dist1, forward2));
		} else if (dist2 < dist1) {
			next_block(curr1, forward1)->split(position(curr1, dist2, forward1));
		} else {
			Leaf* leaf1 = next_block(curr1, forward1);
			Leaf* leaf2 = next_block(curr2, forward2);
			leaf1->join(leaf2, forward1 == forward2);
			
			curr1 = next1;
			curr2 = next2;
		}
	}
	if (not (curr1 == bounds1.second and curr2 == bounds2.second)) {
		print_positions(bounds1, forward1);
		print_positions(bounds2, forward2);
		assert(curr1 == bounds1.second and curr2 == bounds2.second);		
	}
}
