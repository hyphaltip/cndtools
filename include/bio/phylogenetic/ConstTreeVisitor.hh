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

#ifndef __BIO_PHYLOGENETIC_CONSTTREEVISITOR_HH__
#define __BIO_PHYLOGENETIC_CONSTTREEVISITOR_HH__

namespace bio { namespace phylogenetic {

	class Tree;
		
	// A visitor for traversing trees in various orders.  Visitors of
	// this type are NOT allowed to modify the tree during traversal.
	class ConstTreeVisitor {
	public:
		virtual ~ConstTreeVisitor() {}
		virtual void preVisitTree(const Tree* t) {}
		virtual void postVisitTree(const Tree* t) {}
		virtual void preVisitEdge(const Tree* parent, const Tree* child) {}
		virtual void postVisitEdge(const Tree* parent, const Tree* child) {}
	};

} }

#endif // __BIO_PHYLOGENETIC_CONSTTREEVISITOR_HH__
