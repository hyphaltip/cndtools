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

#include "bio/phylogenetic/Tree.hh"
#include "bio/formats/newick.hh"
#include "boost/unordered_set.hpp"
#include "util/options.hh"

using namespace bio;
using namespace formats;
using phylogenetic::Tree;
using boost::unordered_set;

std::string underscoresToSpaces(const std::string& s) {
	std::string str = s;
	std::replace(str.begin(), str.end(), '_', ' ');
	return str;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool keepDescendants = false;
	std::vector<std::string> taxa;

	// Parse command line
	util::options::Parser parser("< newickInput",
								 "Output tree with given taxa as "
								 "contained in newick input");
	parser.addStoreTrueOpt('d', "keep-descendants", 
						   "include in the output tree all descendants of "
						   "given taxa", keepDescendants);
	parser.addAppendArg("taxon", "Taxon to include in the output tree",
						taxa);
	parser.parse(argv, argv + argc);

	try {
		Tree* tree = newick::readTree(std::cin);
		if (tree == NULL) {
			throw std::runtime_error("Invalid input tree");
		}
				
        boost::unordered_set<std::string> taxaSet;
		for (size_t i = 0; i < taxa.size(); ++i) {
			taxaSet.insert(underscoresToSpaces(taxa[i]));
		}
		
		Tree* subtree = tree->getSubtree(taxaSet, keepDescendants);
		
		newick::writeTree(std::cout, subtree);
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
