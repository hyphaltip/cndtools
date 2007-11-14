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

#include "math/LogDouble.hh"
#include "util/stl.hh"
#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;
using math::LogDouble;

#include "OrderedPolygon.hh"
#include "AlignmentReader.hh"

const std::string USAGE = "";

const std::string DESCRIPTION =
	"Given a recombination polygon and a probability distribution over the "
	"parameter space, outputs the posterior probability of each vertex.";

void read_logdoubles(std::istream& stream, std::vector<LogDouble>& lhs) {
	LogDouble ld;
	while (stream >> ld) {
		lhs.push_back(ld);
	}
}

void read_points(std::istream& stream, std::vector< Vector<double> >& points) {
	Vector<double> point(2);
	while (stream >> point) {
		points.push_back(point);
	}
}

Vector<double> fan_ray(Vector<double> point, size_t n) {
	point[0] = std::log(4.0f) + std::log(point[0] / (1.0 - point[0]));
	point[1] = std::log(static_cast<double>(n - 1)) +
		std::log(point[1] / (1.0 - point[1]));
	return point;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	std::string alignment_filename;
	std::string polygon_filename;
	std::string points_filename;
	std::string probs_filename;

	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreArg("alignment_file",
					   "alignment file",
					   alignment_filename);
	parser.addStoreArg("polygon_file",
					   "recombination polygon file",
					   polygon_filename);
	parser.addStoreArg("points_file",
					   "file with points",
					   points_filename);
	parser.addStoreArg("probability_file",
					   "file with probabilities of each point",
					   probs_filename);
	parser.parse(argv, argv + argc);

	try {
		// Read in alignment
		InputFileStream alignment_file(alignment_filename);
		AlignmentReader reader(alignment_file, "", 0);
		size_t n = reader.getAlignment().getNumSeqs() - 1;

		// Read in polygon
		InputFileStream polygon_file(polygon_filename);
		OrderedPolygon<int> p;
		polygon_file >> p;

		// Read in points
		std::vector< Vector<double> > points;
		InputFileStream points_file(points_filename);
		read_points(points_file, points);

		// Read in probabilities
		std::vector<LogDouble> probs;
		InputFileStream probs_file(probs_filename);
		read_logdoubles(probs_file, probs);

		if (probs.size() != points.size()) {
			throw std::runtime_error("Lengths of points and probs files are "
									 "not the same");
		}

		std::vector<LogDouble> posteriors(p.getNumVertices());
		std::vector<size_t> counts(p.getNumVertices());
		LogDouble max_prob;
		size_t max_vnum = 0;
		for (size_t i = 0; i < points.size(); ++i) {
			size_t vnum = p.getVertexNum(fan_ray(points[i], n));
			posteriors[vnum] += probs[i];
			++counts[vnum];
			if (probs[i] > max_prob) {
				max_prob = probs[i];
				max_vnum = vnum;
			}
		}

		LogDouble normalization = util::stl::sum(posteriors.begin(),
												 posteriors.end());

		for (size_t i = 0; i < posteriors.size(); ++i) {
			std::cout << posteriors[i] / normalization << '\t'
					  << counts[i] << '\t'
					  << (i == max_vnum ? 1 : 0) << '\n';
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
