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
#include <list>
#include <vector>
#include <map>
#include <algorithm>

#include "boost/function_output_iterator.hpp"

#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "util/interval.hh"
#include "util/parser.hh"
#include "util/options.hh"

// A coordinate within a Genome
struct GenomicCoord {
	std::string chrom;
	unsigned int coord;

	GenomicCoord(const std::string& chrom, unsigned int coord) 
		: chrom(chrom), coord(coord) {
	}

	bool operator==(const GenomicCoord& x) const {
		return chrom == x.chrom && coord == x.coord;
	}
	bool operator<(const GenomicCoord& x) const {
		return (chrom < x.chrom)
			|| (chrom == x.chrom && coord < x.coord);
	}
	bool operator<=(const GenomicCoord& x) const {
		return (chrom < x.chrom)
			|| (chrom == x.chrom && coord <= x.coord);
	};
};

// A vertex in a graph representing GFF intervals, edges in the graph
// correspond to overlaps between the intervals
struct GFFVertex {
	typedef GenomicCoord coord_type;
	typedef unsigned int difference_type;
	
	bio::gff::GFFRecord rec;
	bool processed;
	bool selected;
	std::list<GFFVertex*> overlaps;

	GFFVertex(const bio::gff::GFFRecord& rec)
		: rec(rec), processed(false), selected(false), overlaps() {
	}

	GenomicCoord start() const {
		return GenomicCoord(rec.getSeqname(), rec.getStart() - 1);
	}
	GenomicCoord end() const {
		return GenomicCoord(rec.getSeqname(), rec.getEnd());
	}
	unsigned long long length() const {
		return rec.getEnd() - rec.getStart() + 1;
	}
	
	bool isProcessed() const { return processed; }
	bool isSelected() const { return selected; }
	void markProcessed() { processed = true; }
	void addOverlap(GFFVertex* v) { overlaps.push_back(v); }
	void select() {
		selected = true;
		std::for_each(overlaps.begin(), overlaps.end(),
					  std::mem_fun(&GFFVertex::markProcessed));
	}
};

// Functor for receiving the results of the overlap analysis of the
// GFF intervals.  It creates an edge in the graph representing each
// overlap that it receives.
struct OverlapRecorder {
 	void operator()(const std::pair<GFFVertex*, GFFVertex*>& overlap) {
		overlap.first->addOverlap(overlap.second);
		overlap.second->addOverlap(overlap.first);
	}
};

// A class for sorting comparisons of GFF vertices.  The base class
// does not do any sorting.
struct Sorter {
	virtual ~Sorter() {}
	virtual bool operator()(const GFFVertex* x, const GFFVertex* y) {
		return false;
	}
};

// Sort GFF vertices by length and then by a secondary sorter
struct LengthSorter : public Sorter {
	Sorter* ss;
	LengthSorter(Sorter* ss) : ss(ss) {}
	bool operator()(const GFFVertex* x, const GFFVertex* y) {
		return x->length() > y->length()
			|| (x->length() == y->length() &&
				(*ss)(x, y));
	}
};

// Sort GFF vertices by start coordinate and then by a secondary sorter
struct StartSorter : public Sorter {
	Sorter* ss;
	StartSorter(Sorter* ss) : ss(ss) {}
	bool operator()(const GFFVertex* x, const GFFVertex* y) {
		return x->start() < y->start()
			|| (x->start() == y->start() && (*ss)(x, y));
	}
};

// A mapping between GFF Source names and priorities.  Larger
// priorities mean higher priority
typedef std::map<std::string, unsigned int> PriorityMap;

// Sort GFF vertices by Source, with higher priority sources coming
// first, and then sort by a secondary sorter
struct SourceSorter : public Sorter {
	Sorter* ss;
	PriorityMap pm;
	
	SourceSorter(Sorter* ss, const PriorityMap& pm) : ss(ss), pm(pm) {}

	bool operator()(const GFFVertex* x, const GFFVertex* y)  {
		return (pm[x->rec.getSource()] > pm[y->rec.getSource()])
			|| (pm[x->rec.getSource()] == pm[y->rec.getSource()] && (*ss)(x, y));
	}
};

// Wraps a pointer to a sorter so that the proper sorting functions
// are called (when the sorter is of a derived sorter class)
struct SorterPtr : public Sorter {
	Sorter* s;
	SorterPtr(Sorter* s) : s(s) {}
	bool operator()(const GFFVertex* x, const GFFVertex* y) {
		return (*s)(x, y);
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Options and defaults
	std::string feature;
	bool verbose = false;
	std::vector<char> sortList;
	std::vector<std::string> sources;
	
	util::options::Parser parser("< gffInput", "");
	parser.addStoreOpt('f', "feature",
					   "feature type for which to detect overlaps",
					   feature);
	parser.addStoreTrueOpt('v', "verbose",
						   "output status information",
						   verbose);
	parser.addAppendConstOpt('l', "length",
							 "give priority to longer features",
							 sortList, 'l');
	parser.addAppendConstOpt('c', "coord",
							 "give priority to features with smaller start "
							 "coordinate",
							 sortList, 'c');
	parser.addAppendConstOpt('s', "source",
							 "give priority to records according to their "
							 "source, as specified by the order of the "
							 "source arguments",
							 sortList, 's');
	parser.addAppendArg("source", "", sources);
	parser.parse(argv, argv + argc);

	// Create the priority map between source names and priorities.
	// All sources listed are given priority > 0.  That way, any
	// source not listed has a default priority of 0.
	PriorityMap pm;
	for (unsigned int i = 0; i < sources.size(); ++i) {
		pm[sources[i]] = sources.size() - i;
	}

	Sorter* nonsorter = new Sorter(); // A sorter that does not sort
	Sorter* sorter = nonsorter; // Start the priority sorter as a nonsorter

	// Compose the sorters backwards to create the proper sort order
	for (std::vector<char>::reverse_iterator it = sortList.rbegin();
		 it != sortList.rend(); ++it) {
		switch (*it) {
		case 'l':
			sorter = new LengthSorter(sorter);
			break;
		case 'c':
			sorter = new StartSorter(sorter);
			break;
		case 's':
			sorter = new SourceSorter(sorter, pm);
			break;
		}
	}

	// Read in the GFF records, creating GFF vertices only for those
	// records that have the specified feature type
	if (verbose) {
		std::cerr << "Reading GFF...\n";
	}
	typedef std::vector<GFFVertex*> VertexList;
	VertexList vertices;
	try {
		// Construct GFF stream for fast reading		
		bio::gff::GFFInputStream gffStream(std::cin);
		bio::gff::GFFRecord rec;
		while (gffStream >> rec) {
			if (feature.empty() || rec.getFeature() == feature) {
				vertices.push_back(new GFFVertex(rec));
			}
		}
	} catch (util::parser::FormatError& err) {
		std::cerr << "Error: " << err.getProblem() << " for line:\n"
				  << err.getLine() << '\n';
		exit(EXIT_FAILURE);
	}

	// Display how many records of the appropriate feature type were
	// found
	if (verbose) {
		std::cerr << "Found "
				  << vertices.size()
				  << (feature.empty() ? "" : " ")
				  << (feature.empty() ? "" : feature)
				  << " records\n";
	}

	// Sort the vertices by start coordinate (necessary precondition
	// for the overlap algorithm)
	if (verbose) {
		std::cerr << "Sorting by start coordinate...\n";
	}
	std::sort(vertices.begin(), vertices.end(), StartSorter(nonsorter));

	// Calculate the overlaps between GFF intervals and make edges in
	// the graph corresponding to the overlaps
	if (verbose) {
		std::cerr << "Calculating overlaps...\n";
	}
	util::interval::overlaps(vertices.begin(), vertices.end(),
							 boost::make_function_output_iterator(OverlapRecorder()));

	// Sort vertices by priority
	if (verbose) {
		std::cerr << "Sorting by specified priorities...\n";
	}
	std::sort(vertices.begin(), vertices.end(), SorterPtr(sorter));

	// Select the set of non-overlapping records by traversing
	// priority-ordered vertex list
	if (verbose) {
		std::cerr << "Selecting set of non-overlapping records...\n";
	}
	VertexList::iterator it;
	for (it = vertices.begin(); it != vertices.end(); ++it) {
		if (!(*it)->isProcessed()) {
			(*it)->select();
			(*it)->markProcessed();
		}
	}

	// Sort vertices by start coordinate before output
	if (verbose) {
		std::cerr << "Sorting by start coordinate...\n";
	}
	std::sort(vertices.begin(), vertices.end(), StartSorter(nonsorter));

	// Write selected records to output
	if (verbose) {
		std::cerr << "Writing GFF...\n";
	}
	unsigned int counter = 0;
	for (it = vertices.begin(); it != vertices.end(); ++it) {
		if ((*it)->isSelected()) {
			std::cout << (*it)->rec;
			++counter;
		}
	}

	// Note how many records were output
	if (verbose) {
		std::cerr << "Wrote " << counter
				  << " non-overlapping"
				  << (feature.empty() ? "" : " ") 
				  << (feature.empty() ? "" : feature)
				  << " records\n";
	}
}
