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

#include "bio/formats/newick.hh"

#include <string>
#include <limits>
#include <stack>

#include "boost/spirit/core.hpp"
#include "boost/spirit/utility/lists.hpp"
#include "boost/spirit/utility/confix.hpp"
using namespace boost::spirit;

#include "bio/phylogenetic/Tree.hh"
#include "bio/phylogenetic/ConstTreeVisitor.hh"
#include "util/io.hh"

namespace bio { namespace formats { namespace newick {

	bool needsQuotes(const std::string& s) {
		return s.find_first_of("\t\n\r\f\v()[]':;,") != std::string::npos;
	}
	
	bool isQuoted(const std::string& s) {
		return s.length() > 1 && s[0] == '\'' && s[s.length() - 1] == '\'';
	}

	std::string quote(const std::string& s) {
		std::string str;
		str.reserve(2 + std::count(s.begin(), s.end(), '\''));
		str += '\'';
		for (std::string::const_iterator it = s.begin(); it != s.end(); ++it) {
			if (*it == '\'') {
				str += '\'';
			}
			str += *it;
		}
		str += '\'';
		return str;
	}
	
	std::string unquote(const std::string& s) {
		std::string us = s;
		std::string::iterator insertPos = us.begin();
		std::string::iterator currPos = us.begin();
		bool useQuote = false;
		while (currPos != us.end()) {
			if (*currPos != '\'' || useQuote) {
				*insertPos = *currPos;
				++insertPos;
				useQuote = false;
			} else {
				useQuote = true;
			}
			++currPos;
		}
		us.resize(insertPos - us.begin());
		return us;
	}

	std::string underscoresToSpaces(const std::string& s) {
		std::string str = s;
		std::replace(str.begin(), str.end(), '_', ' ');
		return str;
	}

	std::string spacesToUnderscores(const std::string& s) {
		std::string str = s;
		std::replace(str.begin(), str.end(), ' ', '_');
		return str;
	}

	
	typedef std::stack<phylogenetic::Tree*> TreeStack;

	struct PushTreeAction {
		TreeStack& s;
		PushTreeAction(TreeStack& s) : s(s) {}
		void operator()(char) const {
			s.push(new phylogenetic::Tree());
		}
		void operator()(char const*, char const*) const {
			s.push(new phylogenetic::Tree());
		}
	};

	struct PopTreeAction {
		TreeStack& s;
		PopTreeAction(TreeStack& s) : s(s) {}
		void operator()(char const*, char const*) const {
			phylogenetic::Tree* descendant = s.top();
			s.pop();
			s.top()->addDescendant(descendant);
		}
	};
	
	struct LabelAction {
		TreeStack& s;
		LabelAction(TreeStack& s) : s(s) {}
		void operator()(char const* first, char const* last) const {
			std::string str(first, last);
			if (isQuoted(str)) {
				str = unquote(str);
			} else {
				str = underscoresToSpaces(str);
			}
			s.top()->setLabel(str);
		}
	};

	struct BranchLengthAction {
		TreeStack& s;
		BranchLengthAction(TreeStack& s) : s(s) {}
		void operator()(const double x) const {
			s.top()->setBranchLength(x);
		}
	};
	
	struct TreeGrammar : public grammar<TreeGrammar>
	{
		TreeStack& stack;
		
		TreeGrammar(TreeStack& stack)
			: stack(stack)
		{}
		
		template <typename ScannerT>
		struct definition
		{
			definition(TreeGrammar const& self)
			{
				tree =
					*comment
					>> descendant_list
					>> *comment
					>> !root_label
					>> *comment
					>> !(':' >> *comment >> branch_length)
					>> *comment
					>> ';';

				descendant_list =
					ch_p('(')[PushTreeAction(self.stack)]
					>> list_p(subtree[PopTreeAction(self.stack)], ',')
					>> ch_p(')');

				subtree =
					*comment
					>> ((descendant_list
						 >> *comment
						 >> !internal_node_label
						 >> *comment
						 >> !(':' >> *comment >> branch_length))
						| (leaf_label
						   >> *comment
						   >> !(':' >> *comment >> branch_length)))
					>> *comment;

				root_label = label[LabelAction(self.stack)];
				internal_node_label = label[LabelAction(self.stack)];
				leaf_label = label[PushTreeAction(self.stack)][LabelAction(self.stack)];
				label = quoted_label | unquoted_label;
				unquoted_label = lexeme_d[+(graph_p - '(' - ')' - '[' - ']' - '\'' - ':' - ';' - ',')];
				quoted_label = lexeme_d[ch_p('\'') >> +((print_p - '\'') | "''") >> ch_p('\'')];
				branch_length = real_p[BranchLengthAction(self.stack)];
				comment = confix_p('[', *anychar_p, ']');
			}

			rule<ScannerT> tree, descendant_list, subtree,
				root_label, internal_node_label, leaf_label,
				label, unquoted_label, quoted_label,
				comment,
				branch_length;
			rule<ScannerT> const&
			start() const { return tree; }
		};

	};

	phylogenetic::Tree* readTree(std::istream& strm) {
		return readTree(util::io::readStream(strm));
	}

	phylogenetic::Tree* readTree(const std::string& str) {
		TreeStack stack;
		TreeGrammar tg(stack);
        parse_info<> info = parse(str.c_str(), tg, space_p);
		return info.hit ? stack.top() : NULL;
	}

	class TreeWriterVisitor : public phylogenetic::ConstTreeVisitor {
	private:
		std::ostream& strm;
		bool addComma;
		bool writeEdgeLengths;
	public:
		TreeWriterVisitor(std::ostream& strm, bool writeEdgeLengths = true)
			: strm(strm), addComma(false), writeEdgeLengths(writeEdgeLengths) {}
		
		void preVisitTree(const phylogenetic::Tree* t) {
			if (!t->isLeaf()) {
				strm << '(';
			}
		}
			
		void postVisitTree(const phylogenetic::Tree* t) {
			if (!t->isLeaf()) {
				strm << ')';
			}
			if (!t->getLabel().empty()) {
				if (needsQuotes(t->getLabel())) {
					strm << quote(t->getLabel());
				} else {
					strm << spacesToUnderscores(t->getLabel());
				}
			}
			if (writeEdgeLengths) {
				strm << ':' << t->getBranchLength();
			}
			addComma = false;
		}
			
		void preVisitEdge(const phylogenetic::Tree* parent,
						  const phylogenetic::Tree* child) {
			if (addComma) {
				strm << ',';
			}
			addComma = false;
		}
			
		void postVisitEdge(const phylogenetic::Tree* parent,
						   const phylogenetic::Tree* child) {
			addComma = true;
		}
	};

	void writeTree(std::ostream& strm,
				   phylogenetic::Tree* t,
				   bool writeEdgeLengths) {
		if (t != NULL) {
			TreeWriterVisitor visitor(strm, writeEdgeLengths);
			t->dfsTraverse(visitor);
		}
		strm << ";\n";
	}

} } }
