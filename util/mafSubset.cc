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

#include "bio/formats/maf.hh"
#include "util/stl.hh"
#include "util/options.hh"
using namespace bio::formats;

struct SequenceSelector {
	SequenceSelector(std::vector<std::string>& keep)
		: genomes(keep.begin(), keep.end()) {}

	bool operator()(const maf::Sequence& s) const {
		return genomes.find(s.getGenomeAndChrom().first) == genomes.end();
	}

	util::stl::hash_set<std::string> genomes;
};

bool is_all_gap_column(std::vector<maf::Sequence>& seqs, size_t col) {
	for (size_t i = 0; i < seqs.size(); ++i) {
		if (seqs[i].text[col] != '-') { return false; }
	}
	return true;
}

void copy_column(std::vector<maf::Sequence>& seqs,
				 size_t source, size_t target) {
	for (size_t i = 0; i < seqs.size(); ++i) {
		seqs[i].text[target] = seqs[i].text[source];
	}
}

void erase_columns(std::vector<maf::Sequence>& seqs,
				   size_t start, size_t end) {
	for (size_t i = 0; i < seqs.size(); ++i) {
		std::string& s = seqs[i].text;
		s.erase(s.begin() + start, s.begin() + end);
	}
}

void remove_gap_columns(std::vector<maf::Sequence>& seqs) {
	if (seqs.empty()) { return; }
	size_t next = 0;
	size_t curr = 0;
	while (curr < seqs.front().text.size()) {
		if (not is_all_gap_column(seqs, curr)) {
			copy_column(seqs, curr, next);
			++next;
		}
		++curr;
	}
	erase_columns(seqs, next, seqs.front().text.size());
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::vector<std::string> genomes;
	
	// Parse command line
	util::options::Parser parser("< mafInput", "Extract subset of MAF "
								 "containing specified genomes");
	parser.addAppendArg("genome", "Name of genome to keep in MAF", genomes);
	parser.parse(argv, argv + argc);

	try {
		maf::InputStream input_stream(std::cin);
		maf::OutputStream output_stream(std::cout, input_stream.getHeader());

		SequenceSelector selector(genomes);
		
		maf::Record rec;
		while (input_stream >> rec) {
			rec.sequences.erase(std::remove_if(rec.sequences.begin(),
											   rec.sequences.end(),
											   selector),
								rec.sequences.end());
			if (rec.sequences.size() >= 2) {
				remove_gap_columns(rec.sequences);
				output_stream << rec;
			}
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
