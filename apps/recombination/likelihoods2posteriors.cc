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
#include "util/string.hh"
#include "filesystem.hh"
using namespace filesystem;
using math::LogDouble;
using util::string::toString;

LogDouble squared_diff(LogDouble x, LogDouble y) {
	if (x < y) { std::swap(x, y); }
	LogDouble diff(x - y);
	return diff * diff;
}

LogDouble stddev(const std::vector<LogDouble>& vals, const LogDouble& mean) {
	LogDouble sum_of_squares;
	for (size_t i = 0; i < vals.size(); ++i) {
		sum_of_squares += squared_diff(vals[i], mean);
	}
	return (sum_of_squares / LogDouble(vals.size() - 1)) ^ 0.5;
}

void read_logdoubles(std::istream& stream, std::vector<LogDouble>& lhs) {
	LogDouble ld;
	while (stream >> ld) {
		lhs.push_back(ld);
	}
}

const std::string USAGE = "";

const std::string DESCRIPTION = "";

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	std::string prior_filename;
	std::vector<std::string> likelihood_filenames;
	
	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreOpt(0, "prior",
					   "File containing prior probability distribution",
					   prior_filename, "FILE");
	parser.addAppendArg("likelihood_file",
						"File containing likelihoods for some data at some "
						"points",
						likelihood_filenames,
						1);
	parser.parse(argv, argv + argc);

	try {
		size_t num_sets = likelihood_filenames.size();

		// Read likelihoods
		std::vector< std::vector<LogDouble> > likelihoods(num_sets);
		size_t n = 0;
		for (size_t i = 0; i < num_sets; ++i) {
			InputFileStream likelihood_file(likelihood_filenames[i]);
			read_logdoubles(likelihood_file, likelihoods[i]);
			if (i == 0) {
				n = likelihoods[i].size();
			} else if (likelihoods[i].size() != n) {
				throw std::runtime_error("Likelihood files not of same length:"
										 + likelihood_filenames[0]
										 + " length = " + toString(n)
										 + " "
										 + likelihood_filenames[i]
										 + "length = "
										 + toString(likelihoods[i].size()));
			}
		}

		// Read prior, if given
		std::vector<LogDouble> prior;
		bool use_prior = not prior_filename.empty();
		if (use_prior) {
			InputFileStream prior_file(prior_filename);
			read_logdoubles(prior_file, prior);
			if (prior.size() != n) {
				throw std::runtime_error("Prior file not of same length as "
										 "likelihood files");
			}
		}
		
		// Calculate summands
		std::vector<LogDouble> summands(n, LogDouble(1));
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < num_sets; ++j) {
				summands[i] *= likelihoods[j][i];
			}
			if (use_prior) {
				summands[i] *= prior[i];
			}
		}

		LogDouble point_sum = util::stl::sum(summands.begin(), summands.end());
		LogDouble mean = point_sum / LogDouble(n);
		LogDouble mean_stddev = stddev(summands, mean) / (LogDouble(n)^0.5);

		std::cerr << "Sample mean: " << mean << '\n'
				  << "Estimate of standard deviation of sample mean: "
				  << mean_stddev << '\n';

		for (size_t i = 0; i < n; ++i) {
			std::cout << summands[i] / mean << '\n';
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
