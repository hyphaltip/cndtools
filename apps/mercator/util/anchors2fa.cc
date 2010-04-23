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
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include "bio/sdb.hh"
#include "bio/formats/fasta.hh"
#include "bio/translation/Table.hh"
#include "util/options.hh"
#include "boost/unordered_map.hpp"

struct Anchor {
	std::string name;
	std::string chrom;
	char strand;
	unsigned long long start;
	unsigned long long end;
	size_t isCoding;
};

inline std::istream& operator>>(std::istream& strm, Anchor& a) {
	return strm >> a.name >> a.chrom >> a.strand >> a.start >> a.end
				>> a.isCoding;
};

struct AnchorSorter {
	bool operator()(const Anchor* a1, const Anchor* a2) const {
		return (a1->chrom < a2->chrom)
			|| (a1->chrom == a2->chrom && a1->start < a2->start);
	}
};

class SequenceExtractor {
public:
	virtual ~SequenceExtractor() {}
	virtual void init(const std::string& filename) = 0;
	virtual std::string getSeq(const std::string& title) = 0;
};

class SDBSequenceExtractor : public SequenceExtractor {
private:
	bio::SDB::DB db;
public:
	void init(const std::string& filename) {
		db.open(filename, true);
	}

	std::string getSeq(const std::string& title) {
		bio::SDB::Record rec;
		db.getRec(title, rec);
		return rec.getSeq();
	}
};

class FASTASequenceExtractor : public SequenceExtractor {
private:
    boost::unordered_map<std::string, std::string> seqMap;
public:
	void init(const std::string& filename) {
		std::ifstream fastaFile(filename.c_str());
		bio::formats::fasta::InputStream fastaStream(fastaFile);
		bio::formats::fasta::Record rec;
		while (fastaStream >> rec) {
			seqMap[rec.title] = rec.sequence;
		}
	}
	
	std::string getSeq(const std::string& title) {
		return seqMap[title];
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Values to be (optionally) specified on the command line
	std::string seqFilename;
	int tableNum = 1;
	SequenceExtractor* seqExtractor = new SDBSequenceExtractor();

	// Parse options
	util::options::Parser parser("< anchorInput",
								 "Output FASTA protein sequence records for "
								 "each genomic anchor interval in the input");
	parser.addStoreOpt('t', "table",
					   "Number of the translation table to use for translation",
					   tableNum, "NUM");
	parser.addStoreConstOpt<SequenceExtractor*>(0, "fasta",
												"Sequence file is a multi-fasta file instead of a "
												"SDB file", seqExtractor,
												new FASTASequenceExtractor());
	parser.addStoreArg("seqFile",
					   "SDB file (or FASTA file if using --fasta option) "
					   "containing genomic sequence corresponding "
					   "to anchors in input",
					   seqFilename);
	parser.parse(argv, argv + argc);

	try {
		seqExtractor->init(seqFilename);
		
		const bio::translation::Table* table =
			bio::translation::Table::getTable(tableNum);
		if (table == NULL) {
			throw std::runtime_error("Invalid translation table: " +
									 util::string::toString(tableNum));
		}

		// Read anchors into a vector
		Anchor a;
		std::vector<Anchor*> anchors;
		while (std::cin >> a) {
			anchors.push_back(new Anchor(a));
		}
		
		// Sort according to coordinate to enable caching of sequence data
		std::sort(anchors.begin(), anchors.end(), AnchorSorter());
		
		// Output sequence for each anchor
		std::string lastChrom = "";
		std::string lastSeq = "";
		bio::formats::fasta::Record faRec;
		bio::formats::fasta::OutputStream fastaOutStream(std::cout);
		std::string dna;
		for (std::vector<Anchor*>::const_iterator pos = anchors.begin();
			 pos != anchors.end(); ++pos) {
			Anchor* a = *pos;
			try {
				if (a->chrom != lastChrom) {
					lastSeq = seqExtractor->getSeq(a->chrom);
					lastChrom = a->chrom;
				}
				if (a->start < 0 ||
					a->end > lastSeq.length() ||
					a->end < a->start) {
					throw std::runtime_error("Invalid coords for anchor "
											 + a->name);
				}
				dna = lastSeq.substr(a->start, a->end - a->start);
				if (a->strand == '-') {
					bio::alphabet::AmbiguousDNA.reverseComplementInPlace(dna);
				}
				
				faRec.title = a->name;
				faRec.sequence = table->translate(dna);
				fastaOutStream << faRec;
			} catch (const std::runtime_error& e) {
				std::cerr << "Error: " << e.what() << '\n';
				continue;
			}
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
