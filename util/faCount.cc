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

#include "bio/formats/fasta/InputStream.hh"
#include "bio/alphabet/AmbiguousNucleotide.hh"
#include "util/string.hh"
#include "util/options.hh"

// This is a C++ rewrite (with my own libraries) of Jim Kent's
// (http://www.cse.ucsc.edu/~kent/) faCount C program

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool showHeader = true;
	bool showTotal = true;

	// Set up options
	util::options::Parser parser("< fastaInput",
								 "Calculates statistics of sequences in fasta records");
	parser.addStoreFalseOpt('h', "no-header", "do not display header",
							showHeader);
	parser.addStoreFalseOpt('t', "no-total", "do not display totals",
							showTotal);
	parser.parse(argv, argv + argc);

	// Integer encodings of nucleotide characters
	int A = bio::alphabet::GenomeDNA.encode('A');
	int C = bio::alphabet::GenomeDNA.encode('C');
	int G = bio::alphabet::GenomeDNA.encode('G');
	int T = bio::alphabet::GenomeDNA.encode('T');
	int N = bio::alphabet::GenomeDNA.encode('N');

	// Initialize counts
	typedef unsigned long long CountType;
	CountType totalLength = 0;
	CountType totalCpgCount = 0;
	std::vector<CountType> baseCounts(bio::alphabet::GenomeDNA.getSize());
	std::vector<CountType> totalBaseCounts(baseCounts.size());

	// Optionally Output header
	if (showHeader) {
		std::cout << "#seq\tlen\tA\tC\tG\tT\tN\tcpg\n";
	}

	// Construct FASTA stream for fast reading
	bio::formats::fasta::InputStream fastaStream(std::cin);
	
	// Step through records
	bio::formats::fasta::Record rec;
    while (fastaStream >> rec) {
		// Reset single record counts to zero
		std::fill(baseCounts.begin(), baseCounts.end(), 0);
        CountType cpgCount = 0;
        int prevBase = -1;

		// Count bases and CpGs
		std::string::const_iterator c;
		for (c = rec.sequence.begin(); c != rec.sequence.end(); ++c) {
			int currBase = bio::alphabet::GenomeDNA.encode(*c);
			++baseCounts[currBase];
			if ((prevBase == C) && (currBase == G)) {
				++cpgCount;
			}
			prevBase = currBase;
		}

		// Output single record counts
		std::cout << rec.title << '\t'
				  << rec.sequence.length() << '\t'
				  << baseCounts[A] << '\t'
				  << baseCounts[C] << '\t'
				  << baseCounts[G] << '\t'
				  << baseCounts[T] << '\t'
				  << baseCounts[N] << '\t'
				  << cpgCount << '\n';

		// Add single record counts to totals
        totalLength += rec.sequence.length();
        totalCpgCount += cpgCount;
        for (unsigned int i = 0; i < baseCounts.size(); i++) {
            totalBaseCounts[i] += baseCounts[i];
		}
	}

	// Optionally show total counts
	if (showTotal) {
		std::cout << "total" << '\t'
				  << totalLength << '\t'
				  << totalBaseCounts[A] << '\t'
				  << totalBaseCounts[C] << '\t'
				  << totalBaseCounts[G] << '\t'
				  << totalBaseCounts[T] << '\t'
				  << totalBaseCounts[N] << '\t'
				  << totalCpgCount << '\n';
	}

	return EXIT_SUCCESS;
}
