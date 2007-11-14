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
#include <string>
#include <vector>

#include "util/stl.hh"
#include "util/string.hh"
#include "util/options.hh"
#include "filesystem.hh"
#include "bio/formats/fasta.hh"
#include "bio/alignment/DNAScoringMatrix.hh"

#include "MaxIdentityAligner.hh"
#include "AnnotatedPolytope.hh"

using bio::alphabet::DNA;
using namespace bio::alignment;
using namespace bio::formats;
using namespace filesystem;
using util::string::toString;

#include <Rational.h>
#include <Matrix.h>
#include <Vector.h>
#include <Set.h>

using polymake::Matrix;
using polymake::Vector;
using polymake::Set;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string match = "0";
	std::string mismatch = "x";
	std::string space = "s";
	std::string gap = "g";
	std::string scoringMatrixFilename;
	std::string seqsFilename;
	size_t start;
	size_t end;

	// Parse command line
	util::options::Parser parser("< fastaInput", "");
	parser.addStoreOpt('m', "match", "match variable", match, "SYMBOL");
	parser.addStoreOpt('x', "mismatch", "mismatch variable", mismatch, "SYMBOL");
	parser.addStoreOpt('s', "space", "space variable", space, "SYMBOL");
	parser.addStoreOpt('g', "gap", "gap variable", gap, "SYMBOL");
	parser.addStoreOpt(0, "matrix", "symbolic scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
	parser.addStoreArg("fastafile", "Sequences for alignment", seqsFilename);
	parser.addStoreArg("start",
					   "Start of interval in which to maximize identities",
					   start);
	parser.addStoreArg("end",
					   "End of interval in which to maximize identities",
					   end);
	parser.parse(argv, argv + argc);

	try {
		// Construct FASTA stream for fast reading
		InputFileStream fastaFile(seqsFilename);
		fasta::InputStream fastaStream(fastaFile);

		// Get first two records
		fasta::Record rec1, rec2;
		fastaStream >> rec1 >> rec2;
		if (not fastaStream) {
			throw std::runtime_error("Alignment requires two FASTA records");
		}

		DNAScoringMatrix<std::string> paramMatrix;
		if (not scoringMatrixFilename.empty()) {
			InputFileStream scoringMatrixFile(scoringMatrixFilename);
			paramMatrix.readMatrix(scoringMatrixFile);
		} else {
			paramMatrix.setMatchScore(match);
			paramMatrix.setMismatchScore(mismatch);
		}

		MaxIdentityAligner<int> aligner(rec1.sequence, rec2.sequence,
										start, end,	
										paramMatrix, space, gap);

		AnnotatedPolytope polytope;
		std::cin >> polytope;

		SetListSection* ftvSection;
		polytope.getSection("FACETS_THRU_VERTICES", ftvSection);
		std::vector< Set<int> > facetsThruVertices = ftvSection->sets;

		MatrixSection<Integer>* rvSection;
		polytope.getSection("REASONABLE_VERTICES", rvSection);
		Matrix<Integer> reasonableVertices = rvSection->matrix;

		MatrixSection<Integer>* fnSection;
		polytope.getSection("FACET_NORMALS", fnSection);
		Matrix<Integer> facetNormals = fnSection->matrix;

		BasicSection* alignmentSection = new BasicSection();
		alignmentSection->title = "CRE_ALIGNMENTS";
		polytope.setSection(alignmentSection);

		MatrixSection<Integer>* identities = new MatrixSection<Integer>();
		identities->title = "CRE_IDENTITIES";
		identities->matrix = Matrix<Integer>(reasonableVertices.rows(), 1);
		polytope.setSection(identities);

		MatrixSection<Integer>* vertexNormals = new MatrixSection<Integer>();
		vertexNormals->title = "VERTEX_NORMALS";
		vertexNormals->matrix = Matrix<Integer>(reasonableVertices.rows(),
												facetNormals.cols());
		polytope.setSection(vertexNormals);
		
		std::vector<std::string> alignments;
		for (int i = 0; i < reasonableVertices.rows(); ++i) {
			const Set<int>& facets = facetsThruVertices[i];
			Vector<Integer> vertexNormal;
			for (Set<int>::const_iterator facetIt = facets.begin();
				 facetIt != facets.end(); ++facetIt) {
				if (facetIt == facets.begin()) {
					vertexNormal = facetNormals[*facetIt];
				} else {
					vertexNormal += facetNormals[*facetIt];
				}
			}
			vertexNormals->matrix.row(i) = vertexNormal;

			std::vector<int> paramVals(vertexNormal.begin() + 1,
									   vertexNormal.end());
			if (paramVals.size() != aligner.getNumParams()) {
				throw std::runtime_error("Dimension mismatch: " +
										 toString(aligner.getNumParams()) +
										 " " +
										 toString(paramVals.size()));
			}

			PairwiseAlignment alignment = aligner.getAlignment(paramVals);
			alignmentSection->lines.push_back(alignment.seq1);
			alignmentSection->lines.push_back(alignment.seq2);
			identities->matrix(i, 0) = static_cast<int>(alignment.getNumIdentities());
		}

		std::cout << polytope;

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
