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
#include <vector>
#include <ctime>

#include "bio/formats/fasta.hh"
#include "util/string.hh"
#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;
using namespace bio::formats;

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int.hpp"
#include "boost/random/variate_generator.hpp"

typedef std::vector<uint32_t> CountVector;
typedef boost::mt19937 base_generator_type;
typedef boost::uniform_int<> distribution_type;
typedef boost::variate_generator<base_generator_type&,
								 distribution_type> variate_generator_type;

base_generator_type generator(std::time(0));
distribution_type uni_dist(0, 3);
variate_generator_type uni(generator, uni_dist);

int encode(char c) {
	c = std::toupper(c);
	switch (c) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		return uni();
	}
}

void count_kmers(const std::string& s, CountVector& v, uint32_t k) {
	if (s.size() < k) { return; }
	const uint32_t mask = (1 << (k * 2)) - 1;
	uint32_t kmer = 0;
	for (size_t i = 0; i < (k - 1); ++i) {
		kmer <<= 2;
		kmer |= encode(s[i]);
	}
	for (size_t i = (k - 1); i < s.size(); ++i) {
		kmer <<= 2;
		kmer &= mask;
		kmer |= encode(s[i]);
		++v[kmer];
	}
}

void count_kmers(std::istream& stream, CountVector& v, uint32_t k) {
	fasta::InputStream fasta_stream(stream);
	fasta::Record rec;
	while (fasta_stream >> rec) {
		count_kmers(rec.sequence, v, k);
	}
}

void mask_seq(std::string& s, size_t start, size_t end) {
	while (start != end) {
		s[start] = std::tolower(s[start]);
		++start;
	}
}

void mask_seq(std::string& s, CountVector& v, uint32_t k, uint32_t d) {
	if (s.size() < k) { return; }
	const uint32_t mask = (1 << (k * 2)) - 1;
	uint32_t kmer = 0;
	for (size_t i = 0; i < (k - 1); ++i) {
		kmer <<= 2;
		kmer |= encode(s[i]);
	}
	for (size_t i = (k - 1); i < s.size(); ++i) {
		kmer <<= 2;
		kmer &= mask;
		kmer |= encode(s[i]);
		if (v[kmer] >= d) {
			mask_seq(s, i - k + 1, i + 1);
		}
	}
}

void mask_stream(std::istream& in_stream, std::ostream& out_stream,
				 CountVector& v, uint32_t k, uint32_t d) {
	fasta::InputStream in_fasta_stream(in_stream);
	fasta::OutputStream out_fasta_stream(out_stream);
	fasta::Record rec;
	while (in_fasta_stream >> rec) {
		mask_seq(rec.sequence, v, k, d);
		out_fasta_stream << rec;
	}
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	size_t k = 14;
	size_t d = 20;
	std::string seq_filename;

	// Parse command line
	util::options::Parser parser("",
								 "Mask k-mers that are overrepresented. "
								 "Outputs a softmasked sequence file.");
	parser.addStoreOpt('k', "", "k-mer size", k);
	parser.addStoreOpt('d', "",
					   "k-mers with duplicity greater than this value "
					   "will be masked", d);
	parser.addStoreArg("fastaFile", "FASTA file containing sequences",
					   seq_filename);
	parser.parse(argv, argv + argc);

	try {
		CountVector counts(1 << (2 * k), 0);
		{
			InputFileStream seq_stream(seq_filename);
			count_kmers(seq_stream, counts, k);
		}
		{
			InputFileStream seq_stream(seq_filename);
			mask_stream(seq_stream, std::cout, counts, k ,d);
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
		
	return EXIT_SUCCESS;
}
