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
#include <fstream>
#include <map>
#include <vector>

#include "bio/formats/agp/Record.hh"
#include "bio/formats/agp/InputStream.hh"
#include "bio/formats/fasta/InputStream.hh"
#include "bio/formats/fasta/OutputStream.hh"
#include "bio/alphabet/AmbiguousNucleotide.hh"
#include "util/stl.hh"
#include "util/options.hh"
#include "filesystem/InputFileStream.hh"
using namespace bio::formats;

struct AGPContigStartSorter {
	bool operator()(const agp::Record* r1,
					const agp::Record* r2) const {
		return r1->getContigStart() < r2->getContigStart();
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string agpFilename;
	
	util::options::Parser parser("< fastaInput > fastaOutput",
								 "Assemble sequences from input into larger "
								 "sequences as specified by an AGP file");
	parser.addStoreArg("apgFile",
					   "AGP file specifying how the sequences in the input "
					   "are assembled into sequences in the output",
					   agpFilename);
	parser.parse(argv, argv + argc);

	// Read input sequences into a map
	typedef util::stl::hash_map<std::string, std::string> SeqMap;
	SeqMap seqRecs;
	bio::formats::fasta::InputStream fastaStream(std::cin);
	bio::formats::fasta::Record fastaRecord;
	while (fastaStream >> fastaRecord) {
		seqRecs[fastaRecord.title] = fastaRecord.sequence;
	}
		
	// Read AGP records into map
	std::map<std::string, std::vector<agp::Record*> > chromRecs;
	filesystem::InputFileStream agpFile(agpFilename);
	//std::ifstream agpFile(agpFilename.c_str());
	agp::InputStream agpStream(agpFile);
	agp::Record rec;
	while (agpStream >> rec) {
		chromRecs[rec.getChrom()].push_back(new agp::Record(rec));
	}

	// Construct each chromosome separately and insert into output
	bio::formats::fasta::OutputStream fastaOutStream(std::cout);
	std::map<std::string, std::vector<agp::Record*> >::iterator it;
	for (it = chromRecs.begin(); it != chromRecs.end(); ++it) {
		const std::string& chrom = it->first;
		std::vector<agp::Record*>& recs = it->second;

		// Order the AGP records by contig start coordinate
		std::sort(recs.begin(), recs.end(), AGPContigStartSorter());

		// Determine the chrom length and make a new sequence
		std::string chromSeq(recs.back()->getContigEnd(), 'N');
		std::vector<agp::Record*>::iterator recIt;
		for (recIt = recs.begin(); recIt != recs.end(); ++recIt) {
			if ((*recIt)->isGap()) {
				continue;
			}

			std::string sourceName = (*recIt)->getSourceAccession();
			SeqMap::iterator seqIt = seqRecs.find(sourceName);
			if (seqIt == seqRecs.end()) {
				std::cerr << "Error: "
						  << "Source record not found in input: "
						  << sourceName << '\n';
				exit(EXIT_FAILURE);
			}
			std::string seq =
				seqIt->second.substr((*recIt)->getSourceStart() - 1,
									 (*recIt)->getSourceEnd() -
									 (*recIt)->getSourceStart() + 1);
			if ((*recIt)->getSourceOrientation() == '-') {
				bio::alphabet::AmbiguousDNA.reverseComplementInPlace(seq);
			}

			chromSeq.replace((*recIt)->getContigStart() - 1, seq.size(), seq);
		}

		fastaOutStream << bio::formats::fasta::Record(chrom, chromSeq);
	}
	
	return EXIT_SUCCESS;
}
