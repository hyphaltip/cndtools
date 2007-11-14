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

#include "polytope/Polytope.hh"
#include "polytope/formats/polymake/InputStream.hh"
#include "bio/alignment/AffineGapNWPairwiseConsensusAligner.hh"
#include "bio/alignment/AmbiguousDNAScoringMatrix.hh"
#include "bio/formats/fasta.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "filesystem.hh"
using namespace bio::alignment;
using namespace filesystem;
using namespace bio::formats;
using namespace polytope;
using util::string::toString;

typedef int T;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	std::string polytopeFilename;
	
	// Set up option parser
	util::options::Parser parser("> rays",
								 "Outputs for each vertex of the input "
								 "polytope a parameter ray for which that "
								 "vertex is optimal");
	
	parser.addStoreArg("polytopeFile",
					   "POLYMAKE polytope file",
					   polytopeFilename);
	parser.parse(argv, argv + argc);

	try {
		InputFileStream polytopeFile(polytopeFilename);
		formats::polymake::InputStream polymakeInputStream(polytopeFile);
		Polytope<T> p;
		polymakeInputStream >> p;
		
		const Polytope<T>::FacetNormalList facetNormals = p.getFacetNormals();
		const Polytope<T>::FacetsThruVertexList facetsThruVertices = p.getFacetsThruVertices();

		for (size_t i = 0; i < p.getNumVertices(); ++i) {
			Polytope<T>::Ray paramRay(p.getAmbientDim() + 1);
			typedef Polytope<T>::FacetsThruVertex::const_iterator FacetIndexIterator;
			for (FacetIndexIterator fi = facetsThruVertices[i].begin();
				 fi != facetsThruVertices[i].end(); ++fi) {
				paramRay += facetNormals[*fi];
			}
			std::cout << paramRay << '\n';
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
