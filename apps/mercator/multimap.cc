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

#include "multimap.hh"
#include "anchor.hh"
#include "chromosome.hh"
#include "assembled.hh"
#include "edge.hh"
#include "genome.hh"
#include "clique.hh"
#include "run.hh"
#include "mask.hh"

#include "util/string.hh"
#include "filesystem/Path.hh"
#include "util/stl.hh"

void unflipAnchors() {
	for (size_t i = 0; i < Genome::getNumGenomes(); ++i) {
		Genome::getGenome(i)->unflipAnchors();
	}
}

void checkCliques() {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		Genome* gen = Genome::getGenome(g);
		for (size_t c = 0; c < gen->getNumChroms(); ++c) {
			Clique* curr = gen->getChrom(c)->getFirstClique();
			while (curr != NULL) {
				if (!curr->isConnected()) {
					std::cerr << "Found unconnected clique";
					assert(false);
				}
				curr = curr->nextClique(g);
			}
		}
	}
}

void removeAllCliques() {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		Genome* gen = Genome::getGenome(g);
		for (size_t c = 0; c < gen->getNumChroms(); ++c) {
			Clique* curr = gen->getChrom(c)->getFirstClique();
			while (curr != NULL) {
				Clique* next = curr->nextClique(g);
				delete curr;
				curr = next;
			}
		}
	}
}

void removeRuns() {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		Genome* gen = Genome::getGenome(g);
		for (size_t c = 0; c < gen->getNumChroms(); ++c) {
			for (Clique* curr = gen->getChrom(c)->getFirstClique();
				 curr != NULL;
				 curr = curr->nextClique(g)) {				 
				if (curr->isInRun()) {
					delete curr->getRun();
				}
			}
		}
	}
}

int countCliquesInRuns(vector<Run*>& runs) {
	int numCliques = 0;
	vector<Run*>::const_iterator run;
	for (run = runs.begin(); run != runs.end(); ++run) {
		numCliques += (*run)->getLength();
	}
	return numCliques;
}

size_t countRepetitiveAnchors() {
	const vector<Genome*>& genomes = Genome::getGenomes();
	return util::stl::sum(genomes.begin(), genomes.end(),
						  std::mem_fun(&Genome::getNumAnchorsRepetitive));
}

void printRuns(const vector<Run*>& runs,
			   ostream& strm) {
	vector<Run*>::const_iterator run;
	for (run = runs.begin(); run != runs.end(); ++run) {
		(*run)->dump(strm);
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if (g) { strm << '\t'; }
			strm << "NA";
		}
		strm << '\n';
	}
}

void debugRuns() {
	vector<Run*> runs;
	getRuns(runs);
	cerr << "Runs:\n";
	for (size_t i = 0; i < runs.size(); ++i) {
		cerr << *runs[i];
		cerr << '\n';
	}
}

void printMap(const vector<Run*>& runs, ostream& strm, const bool extend) {
	for (size_t i = 0; i < runs.size(); ++i) {
		runs[i]->printMapLine(strm, extend);
	}
}

void printPairwiseHits(const vector<Run*>& runs, ostream& strm) {
	for (size_t i = 0; i < runs.size(); ++i) {
		runs[i]->printPairwiseHits(strm);
	}
}

void markRepeats(const size_t repeatNum,
				 const float repeatPct) {
	for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
		Genome::getGenome(g)->markRepeats(repeatNum, repeatPct);
	}
}

Clique* findClique(Anchor* startAnchor,
				   const Mask& notAllowed,
				   const bool incomplete=false) {
	if (startAnchor->isInRun() || startAnchor->isMarked()) {
		return NULL;
	}

	Clique c;
	c.addAnchor(startAnchor);
	
	// Mask to indicate which genomes we have looked at for this
	// clique
	Mask examined(0);

	while (c.getMask() != examined) {
		size_t f = firstInMask(c.getMask() ^ examined);
					
		Anchor* a = c.getAnchor(f);
		assert(!a->isInRun());
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if (g == f) {
				continue;
			}

			if (a->hasEdgesTo(g)) {
				if (notAllowed.test(g)) {
					return NULL;
				}
				
				Anchor* bestHit = a->getBestEdge(g)->getOtherAnchor(a);

				if (bestHit->isMarked()) {
					return NULL;
				} else if (c.hasGenome(g) && c.getAnchor(g) == bestHit) {
					continue;
				} else if (!c.hasGenome(g) &&
						   (incomplete || a == startAnchor)) {
					c.addAnchor(bestHit);
				} else {
					return NULL;
				}
			} else if (c.hasGenome(g) && !incomplete) {
				return NULL;
			}
		}
		examined.set(f);
	}

	return new Clique(c);
}

size_t findCliques(const size_t minSize,
				   bool incomplete=false) {
	size_t numCliquesAdded = 0;
	Mask processed(0);

	for (size_t base = 0; base < Genome::getNumGenomes(); ++base) {
		Genome* baseGen = Genome::getGenome(base);
		
		for (size_t chrom = 0; chrom < baseGen->getNumChroms(); ++chrom) {
			for (Anchor* curr = baseGen->getChrom(chrom)->getFirstAnchor();
				 curr != NULL;
				 curr = curr->nextAnchor()) {

				Clique* c = findClique(curr, processed, incomplete);

				if (c == NULL) {
					continue;
				} else if (c->getSize() < minSize) {
					delete c;
				} else {
					++numCliquesAdded;
					c->useClique();
					Run* r = new Run();
					r->addCliqueToRight(c);
					r->claimCliques();
				}
			}
		}
		
		processed.set(base);
	}
	
	return numCliquesAdded;
}

void findCliquesIter(const size_t minSize,
					 const bool incomplete=false) {
	
	size_t numCliques = 0;
	
	// Repeat until we stop adding cliques
	size_t iter = 1;
	while (true) {
		size_t added = findCliques(minSize, incomplete);

		if (added == 0) {
			break;
		}

		numCliques += added;

		cerr << "Iteration " << iter << ": "
			 << "Added " << added
			 << " cliques of size >=" << minSize
			 << " (Total cliques: " << numCliques << ")" << '\n';
		
		++iter;
	}
}

RunStats calcRunStats(vector<Run*>& runs) {
	RunStats stats;

	vector<Run*> myruns(runs);
	sort(myruns.begin(), myruns.end(), RunLengthSorter());

	stats.numRuns = myruns.size();
	
	// Calculate the number of cliques
	stats.numCliques = 0;
	vector<Run*>::const_iterator r;
	for (r = myruns.begin(); r != myruns.end(); ++r) {
		stats.numCliques += (*r)->getLength();
	}

	if (myruns.empty()) {
		stats.meanRunLength = 0;
		stats.medianRunLength = 0;
		stats.maxRunLength = 0;
		stats.minRunLength = 0;
	} else {
		if (myruns.size() % 2 == 0) {
			stats.medianRunLength = ((myruns[stats.numRuns / 2]->getLength() +
									  myruns[stats.numRuns / 2 - 1]->getLength()) /
									 2);
		} else {
			stats.medianRunLength = myruns[stats.numRuns / 2]->getLength();
		}
		
		stats.maxRunLength = myruns.back()->getLength();
		stats.minRunLength = myruns.front()->getLength();
		stats.meanRunLength = (static_cast<double>(stats.numCliques)
							   / stats.numRuns);
	}

	return stats;
}

void printCounts() {
	vector<Run*> runs;
	getRuns(runs);
	// Calculate number of cliques in runs
	cerr << "Number of runs: " << runs.size()
		 << " (using " << countCliquesInRuns(runs) << " cliques)" << '\n';
}

void dumpHits(const size_t num,
			  const filesystem::Path& outDir) {
	std::ofstream hitFile;
	(outDir / (util::string::toString(num) + ".hits")).openForOutput(hitFile);
	Genome::writeActiveHits(hitFile, 0, 1);
}

void dumpRuns(const size_t num,
			  const filesystem::Path& outDir) {
	std::ofstream runFile;
	(outDir / (util::string::toString(num) + ".runs")).openForOutput(runFile);
	vector<Run*> runs;
	getRuns(runs);
	orderRuns(runs);
	printRuns(runs, runFile);
}
			  
// Strategy:
// 1. Mark repetitive anchors (parameters?)
// 2. Find best cliques (max E-value parameter?, iterate?)
// 3. Join (maximum distance parameter)
// 4. Filter edges within significant runs

void makeMap(vector<Run*>& runs, const MapOptions& options) {
	size_t hitnum = 0;
	size_t runnum = 0;

	// Remove cliques and runs from possible previous run
	removeRuns();
	removeAllCliques();

	// Initialize edges
	cerr << "Initializing edges...\n";
	Genome::initEdges(options.prunePct);

	if (options.outputHits) { dumpHits(hitnum++, options.outdir); }

	// Set minimum run length
	Run::setMinRunLength(options.minRunLength);

	// Set padding
	Assembled::setPadding(options.padding);
	
	for (size_t size = Genome::getNumGenomes(); size > 1; --size) {
		cerr << "Marking repetitive anchors...\n";
		markRepeats(options.repeatNum, options.repeatPct);
		std::cerr << "Number of repetitive anchors: "
				  << countRepetitiveAnchors() << '\n';
		
		cerr << "Finding cliques of minimum size " << size << '\n';
		findCliquesIter(size);
		printCounts();
		
		cerr << "Joining runs with maximum distance " << options.maxDist << '\n';
		joinRuns(size, options.maxDist);
		printCounts();
		if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
		
		cerr << "Filtering intrarun edges...\n";
		filterIntraRunEdges();
		breakRuns();

		cerr << "Checking cliques...\n";
		checkCliques();
		
		if (options.outputHits) { dumpHits(hitnum++, options.outdir); }
	}

	cerr << "Marking repetitive anchors...\n";
	markRepeats(options.repeatNum, options.repeatPct);
	std::cerr << "Number of repetitive anchors: "
			  << countRepetitiveAnchors() << '\n';
	
	cerr << "Finding cliques of minimum size " << 2 << '\n';
	findCliquesIter(2, true);
	printCounts();
	
	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }

	cerr << "Checking cliques...\n";
	checkCliques();
	
	if (options.minRunLength > 1) {
		cerr << "Removing singletons...\n";
		removeSingletons();
		printCounts();

		cerr << "Checking cliques...\n";
		checkCliques();
		
		cerr << "Joining runs without maximum distance...\n";
		joinRuns();
		printCounts();
		if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
	}

	cerr << "Checking cliques...\n";
	checkCliques();
	
	cerr << "Removing non-significant runs...\n";
	removeInsignificantRuns();
	printCounts();

	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
	
	cerr << "Filtering intrarun edges...\n";
	filterIntraRunEdges();
	cerr << "Filtering interrun edges...\n";
	filterInterRunEdges();

	if (options.outputHits) { dumpHits(hitnum++, options.outdir); }
	breakRuns();

	cerr << "Checking cliques...\n";
	checkCliques();
	
	cerr << "Marking repetitive anchors...\n";
	markRepeats(options.repeatNum, options.repeatPct);
	std::cerr << "Number of repetitive anchors: "
			  << countRepetitiveAnchors() << '\n';

	unflipAnchors();
	
	cerr << "Finding cliques of minimum size " << 2 << '\n';
	findCliquesIter(2, true);
	printCounts();

	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }	
	
	cerr << "Removing non-significant runs...\n";
	removeInsignificantRuns();
	printCounts();

	if (options.minRunLength > 1) {
		cerr << "Removing singletons...\n";
		removeSingletons();
		printCounts();
	}

	cerr << "Checking cliques...\n";
	checkCliques();
	
	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
	
	// Break runs that create cycles in draft genomes
	cerr << "Breaking cyclic runs...\n";
	breakCyclicRuns();
	printCounts();

	// Assemble draft genomes
	cerr << "Assembling draft genomes...\n";
	Genome::assembleDraftGenomes();
	printCounts();

	cerr << "Checking cliques...\n";
	checkCliques();
	
	// Collect final runs
	getRuns(runs);
	
	// Order runs
	orderRuns(runs);

	// Number runs
	numberRuns(runs);
	
	cerr << "Map-making completed\n";
}

void joinPairwiseMaps(vector<Run*>& runs, const MapOptions& options) {
	size_t hitnum = 0;
	size_t runnum = 0;

	// Remove cliques and runs from possible previous run
	removeRuns();
	removeAllCliques();

	// Initialize edges
	cerr << "Initializing edges...\n";
	Genome::initEdges(options.prunePct);

	if (options.outputHits) { dumpHits(hitnum++, options.outdir); }

	// Set minimum run length
	Run::setMinRunLength(options.minRunLength);

	// Set padding
	Assembled::setPadding(options.padding);


	cerr << "Marking repetitive anchors...\n";
	markRepeats(options.repeatNum, options.repeatPct);
	std::cerr << "Number of repetitive anchors: "
			  << countRepetitiveAnchors() << '\n';
	
	cerr << "Finding cliques of minimum size " << 2 << '\n';
	findCliquesIter(2);
	printCounts();
	
	cerr << "Joining runs with maximum distance " << options.maxDist << '\n';
	joinRuns(1, options.maxDist);
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
	
	cerr << "Filtering intrarun edges...\n";
	filterIntraRunEdges();
	breakRuns();
	
	cerr << "Checking cliques...\n";
	checkCliques();
	
	if (options.outputHits) { dumpHits(hitnum++, options.outdir); }

	cerr << "Marking repetitive anchors...\n";
	markRepeats(options.repeatNum, options.repeatPct);
	std::cerr << "Number of repetitive anchors: "
			  << countRepetitiveAnchors() << '\n';
	
	cerr << "Finding cliques of minimum size " << 2 << '\n';
	findCliquesIter(2, true);
	printCounts();
	
	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }	

	cerr << "Checking cliques...\n";
	checkCliques();
	
	if (options.minRunLength > 1) {
		cerr << "Removing singletons...\n";
		removeSingletons();
		printCounts();

		cerr << "Checking cliques...\n";
		checkCliques();
		
		cerr << "Joining runs without maximum distance...\n";
		joinRuns();
		printCounts();
		if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }		
	}

	cerr << "Checking cliques...\n";
	checkCliques();
	
	cerr << "Removing non-significant runs...\n";
	removeInsignificantRuns();
	printCounts();

	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
	
	cerr << "Filtering intrarun edges...\n";
	filterIntraRunEdges();
	cerr << "Filtering interrun edges...\n";
	filterInterRunEdges();

	if (options.outputHits) { dumpHits(hitnum++, options.outdir); }
	breakRuns();

	cerr << "Checking cliques...\n";
	checkCliques();
	
	cerr << "Marking repetitive anchors...\n";
	markRepeats(options.repeatNum, options.repeatPct);
	std::cerr << "Number of repetitive anchors: "
			  << countRepetitiveAnchors() << '\n';
	
	cerr << "Finding cliques of minimum size " << 2 << '\n';
	findCliquesIter(2, true);
	printCounts();

	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }	
	
	cerr << "Removing non-significant runs...\n";
	removeInsignificantRuns();
	printCounts();

	if (options.minRunLength > 1) {
		cerr << "Removing singletons...\n";
		removeSingletons();
		printCounts();
	}

	cerr << "Checking cliques...\n";
	checkCliques();
	
	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
	
	// Break runs that create cycles in draft genomes
	cerr << "Breaking cyclic runs...\n";
	breakCyclicRuns();
	printCounts();

	// Assemble draft genomes
	cerr << "Assembling draft genomes...\n";
	Genome::assembleDraftGenomes();
	printCounts();

	cerr << "Checking cliques...\n";
	checkCliques();
	
	// Collect final runs
	getRuns(runs);
	
	// Order runs
	orderRuns(runs);

	// Number runs
	numberRuns(runs);
	
	cerr << "Map-making completed\n";
}

void joinPairwiseMaps2(vector<Run*>& runs, const MapOptions& options) {
	size_t hitnum = 0;
	size_t runnum = 0;

	// Remove cliques and runs from possible previous run
	removeRuns();
	removeAllCliques();

	// Initialize edges
	cerr << "Initializing edges...\n";
	Genome::initEdges(options.prunePct);

	if (options.outputHits) { dumpHits(hitnum++, options.outdir); }

	// Set minimum run length
	Run::setMinRunLength(options.minRunLength);

	// Set padding
	Assembled::setPadding(options.padding);
	
	cerr << "Finding cliques of minimum size " << 2 << '\n';
	findCliquesIter(2, true);
	printCounts();
		
	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }	

	cerr << "Checking cliques...\n";
	checkCliques();
	
	if (options.minRunLength > 1) {
		cerr << "Removing singletons...\n";
		removeSingletons();
		printCounts();

		cerr << "Checking cliques...\n";
		checkCliques();
		
		cerr << "Joining runs without maximum distance...\n";
		joinRuns();
		printCounts();
		if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }		
	}

	cerr << "Checking cliques...\n";
	checkCliques();
	
	cerr << "Removing non-significant runs...\n";
	removeInsignificantRuns();
	printCounts();

	cerr << "Joining runs without maximum distance...\n";
	joinRuns();
	printCounts();
	if (options.outputRuns) { dumpRuns(runnum++, options.outdir); }
	
	// Break runs that create cycles in draft genomes
	cerr << "Breaking cyclic runs...\n";
	breakCyclicRuns();
	printCounts();

	// Assemble draft genomes
	cerr << "Assembling draft genomes...\n";
	Genome::assembleDraftGenomes();
	printCounts();

	cerr << "Checking cliques...\n";
	checkCliques();
	
	// Collect final runs
	getRuns(runs);
	
	// Order runs
	orderRuns(runs);

	// Number runs
	numberRuns(runs);
	
	cerr << "Map-making completed\n";
}
