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
#include <sstream>
#include <vector>
#include <map>

#include "bio/formats/fasta.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "filesystem.hh"
using namespace filesystem;

#include "SegmentTreeSet.hh"

typedef SegmentTreeSet::Leaf Segment;
typedef std::vector<Segment*> SegList;
typedef std::map<std::string, SegList> GenomeSegs;
typedef util::stl::hash_map<const Segment*, size_t> SegmentLabels;
typedef std::pair<std::string, std::string> StringPair;
typedef std::map<StringPair, OutputFileStream*> HitFileMap;

// Separates the genome and chromosome parts of sequence identifier
StringPair genome_and_chrom(const std::string& s) {
	std::string::size_type dot_pos = s.find('.');
	if (dot_pos == std::string::npos) {
		throw std::runtime_error("Sequence name not of the form GENOME.CHROM: "
								 + s);
	}
	return StringPair(s.substr(0, dot_pos), s.substr(dot_pos + 1));
}

// Sorts segments by interval order
struct SegSorter {
	bool operator()(const Segment* s1, const Segment* s2) const {
		return s1->getInterval() < s2->getInterval();
	}
};

// Group segments by genome
void make_genome_segs(GenomeSegs& genome_segs,
					  SegList& segs,
					  size_t min_length) {
	typedef SegList::const_iterator SegListIter;
	for (SegListIter s = segs.begin(); s != segs.end(); ++s) {
		Segment* seg = *s;
		if (not seg->isUnique() and seg->getLength() >= min_length) {
			StringPair names = genome_and_chrom(seg->getSeqName());
			genome_segs[names.first].push_back(seg);
		}
	}
}

// Sort segments within each genome and assign labels in increasing order
void make_seg_labels(SegmentLabels& labels, GenomeSegs& genome_segs) {
	typedef GenomeSegs::iterator GenomeSegsIter;
	typedef SegList::const_iterator SegListIter;
	for (GenomeSegsIter g = genome_segs.begin(); g != genome_segs.end(); ++g) {
		SegList& seg_list = g->second;
		std::sort(seg_list.begin(), seg_list.end(), SegSorter());
		size_t label = 0;
		for (SegListIter s = seg_list.begin(); s != seg_list.end(); ++s) {
			labels[*s] = label++;
		}
	}
}

// Write anchor segment intervals to files for each genome
void write_anchors(const GenomeSegs& genome_segs,
				   const SegmentLabels& labels,
				   const Path& out_dir,
				   const std::string anchor_ext) {
	typedef GenomeSegs::const_iterator GenomeSegsIter;
	typedef SegList::const_iterator SegListIter;
	for (GenomeSegsIter g = genome_segs.begin(); g != genome_segs.end(); ++g) {
		const std::string genome = g->first;
		const SegList& seg_list = g->second;
		OutputFileStream anchor_file(out_dir / (genome + anchor_ext));
		for (SegListIter s = seg_list.begin(); s != seg_list.end(); ++s) {
			bio::genome::BasicInterval i((*s)->getInterval());
			i.setChrom(genome_and_chrom(i.getChrom()).second);
			anchor_file << labels.find(*s)->second << '\t'
						<< i.getChrom() << '\t'
						<< i.getStrand() << '\t'
						<< i.getStart() << '\t'
						<< i.getEnd() << '\t'
						<< 0 << '\n'; // non-coding
		}
	}
}

// Create hit files for all pairs of genomes
void open_hit_files(const GenomeSegs& genome_segs,
					const Path& out_dir,
					const std::string hit_ext,
					HitFileMap& hit_files,
					bool output_self_hits) {
	typedef GenomeSegs::const_iterator Iter;
	for (Iter g1 = genome_segs.begin(); g1 != genome_segs.end(); ++g1) {
		for (Iter g2 = genome_segs.begin(); g2 != genome_segs.end(); ++g2) {
			if (g1->first > g2->first
				or (not output_self_hits and g1->first == g2->first)) {
				continue;
			}
			StringPair p(g1->first, g2->first);
			Path hit_filename = out_dir / (p.first + "-" + p.second + hit_ext);
			hit_files[p] = new OutputFileStream(hit_filename);
		}
	}
}

// Close all open hit files
void close_hit_files(HitFileMap& hit_files) {
	typedef HitFileMap::iterator HitFileIt;
	for (HitFileIt it = hit_files.begin(); it != hit_files.end(); ++it) {
		it->second->close();
	}
}

// Write a hit between anchor S1 and anchor S2
void write_hit(const Segment* s1,
			   const Segment* s2,
			   const SegmentLabels& labels,
			   HitFileMap& hit_files,
			   bool output_self_hits) {
	std::string genome1 = genome_and_chrom(s1->getSeqName()).first;
	std::string genome2 = genome_and_chrom(s2->getSeqName()).first;
	if (genome1 > genome2) {
		std::swap(s1, s2);
		std::swap(genome1, genome2);
	} else if (not output_self_hits and genome1 == genome2) {
		return;
	}
	std::ostream& hit_file = *hit_files[StringPair(genome1, genome2)];
	hit_file << labels.find(s1)->second << '\t'
			 << labels.find(s2)->second << '\t'
			 << s1->getLength() << '\t'
			 << 0 << '\n';
}

void orient_trees(SegmentTreeSet::TreeSet& trees) {
	typedef SegmentTreeSet::TreeSet::iterator TreeIter;
	for (TreeIter t = trees.begin(); t != trees.end(); ++t) {
		(*t)->orientTree();
	}
}

// Write all hits represented by the segment trees
void write_hits(const SegmentTreeSet::TreeSet& trees,
				const SegmentLabels& labels,
				size_t min_length,
				HitFileMap& hit_files,
				bool output_self_hits) {
	typedef SegmentTreeSet::TreeSet::const_iterator TreeIter;
	for (TreeIter t = trees.begin(); t != trees.end(); ++t) {
		if ((*t)->getLength() < min_length) { continue; }
		SegList segments;
		(*t)->getLeaves(segments);
		for (size_t i = 0; i < segments.size(); ++i) {
			for (size_t j = 0; j < i; ++j) {
				write_hit(segments[i], segments[j], labels, hit_files,
						  output_self_hits);
			}
		}
	}
}

// Add sequences (and lengths) to the segment tree set
void add_seqs(SegmentTreeSet& sts, std::istream& stream) {
	std::string name;
	SegmentTreeSet::Position length;
	while (stream >> name >> length) {
		sts.addSeq(name, length);
	}
}

// Add hits to the segment tree set from MUMmer output
void add_hits(SegmentTreeSet& sts,
			  std::istream& stream,
			  bool output_self_hits) {
	std::string ref, query, ref_strand, query_strand;
	size_t ref_start, query_start, length;

	while (stream >> ref >> ref_strand >> ref_start
		   >> query >> query_strand >> query_start >> length) {
		// Check that we have seen these sequence names before
		if (not sts.hasSeq(ref)) {
			throw std::runtime_error("Sequence in match file not listed in "
									 "sequence length file: " + ref);
		}
		if (not sts.hasSeq(query)) {
			throw std::runtime_error("Sequence in match file not listed in "
									 "sequence length file: " + query);
		}

		// Check for hits within the same genome
		std::string ref_genome = genome_and_chrom(ref).first;
		std::string query_genome = genome_and_chrom(query).first;
		if (ref_genome == query_genome and not output_self_hits) {
			continue;
		}

		sts.addMatch(ref, ref_start - 1, ref_strand == "+",
					 query, query_start - 1, query_strand == "+",
					 length);
	}
}

const std::string DESCRIPTION = "";
const std::string USAGE = "";

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string seq_length_filename;
	std::string match_filename;
	size_t min_length = 20;
	std::string out_dir = ".";
	std::string anchor_ext = ".mum.anchors";
	std::string hit_ext = ".mum.hits";
	bool output_self_hits = false;
	
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreOpt('l', "min-length", "Minimum anchor length",
					   min_length, "INT");
	parser.addStoreOpt('o', "out-dir", "Output directory", out_dir, "DIR");
	parser.addStoreOpt(0, "anchor-ext", "Extension for anchor files",
					   anchor_ext, "STRING");
	parser.addStoreOpt(0, "hit-ext", "Extension for hit files",
					   hit_ext, "STRING");
	parser.addStoreTrueOpt(0, "self-hits", "Output hits within a genome",
						   output_self_hits);
	parser.addStoreArg("seqLengthFile",
					   "File containing sequence names and lengths",
					   seq_length_filename);
	parser.addStoreArg("matchFile",
					   "File containing matches between sequences",
					   match_filename);
	parser.parse(argv, argv + argc);

	try {
		InputFileStream seq_length_file(seq_length_filename);
		InputFileStream match_file(match_filename);

		SegmentTreeSet sts;
		SegmentTreeSet::TreeSet trees;
		SegList all_segs;
		GenomeSegs genome_segs;
		SegmentLabels labels;
		HitFileMap hit_files;

		add_seqs(sts, seq_length_file);
		add_hits(sts, match_file, output_self_hits);
		sts.getTrees(trees);
		sts.getSegments(all_segs);
		make_genome_segs(genome_segs, all_segs, min_length);
		make_seg_labels(labels, genome_segs);
		orient_trees(trees);
		write_anchors(genome_segs, labels, out_dir, anchor_ext);
		open_hit_files(genome_segs, out_dir, hit_ext, hit_files,
					   output_self_hits);
		write_hits(trees, labels, min_length, hit_files, output_self_hits);
		close_hit_files(hit_files);
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
