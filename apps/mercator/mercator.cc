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

#include <ctime>

#include "util/string.hh"
#include "util/options.hh"

#include "types.hh"
#include "multimap.hh"
#include "genome.hh"

void writeVersion(std::ostream& strm) {
	strm << "mapmaker " << BUILD_VERSION << " (" << BUILD_DATE << ")" << '\n';
}

int main(int argc, const char* argv[]) {
	// Values to be (optionally) specified on the command line
	MapOptions options;
	options.repeatNum = 2;
	options.repeatPct = 0.90;
	options.maxE = 1;
	options.prunePct = 0.8;
	options.maxDist = 300000;
	options.minRunLength = 2;
	options.padding = 100;
	options.indir = ".";
	options.outdir = ".";
	options.quiet = false;
	options.outputHits = false;
	options.outputRuns = false;
	options.phitFilename = "";
	
	std::vector<std::string> nondraftGenomes;
	std::vector<std::string> draftGenomes;

	bool version = false;

	util::options::Parser parser("[[-d|--draft] GENOME]*", "");

	parser.addStoreOpt(0, "repeat-num",
					   "delay mapping of anchors with at least NUM high-scoring hits to any other genome",
					   options.repeatNum, "NUM");
	parser.addStoreOpt(0, "repeat-pct",
					   "used in conjunction with the \"repeat-num\" option.  \"high-scoring\" hits are considered to be those with scores at least PCT of the maximum hit score", 
					   options.repeatPct, "PCT");
	parser.addStoreOpt('e', "max-eval",
					   "filter out hits with e-value more than EVAL",
					   options.maxE, "EVAL");
	parser.addStoreOpt(0, "prune-pct",
					   "hits with scores less than PCT of the maximum hit score are filtered out",
					   options.prunePct, "PCT");
	parser.addStoreOpt('j', "join-distance", 
					   "max distance at which to join cliques",
					   options.maxDist, "DISTANCE");
	parser.addStoreOpt('l', "min-run-length",
					   "minimum number of consistent cliques for a run",
					   options.minRunLength, "NUM");
	parser.addStoreOpt(0, "padding",
					   "number of Ns to put in between assembled contigs",
					   options.padding, "NUM");
	parser.addAppendOpt('d', "draft",
						"name of a genome to be considered as a draft genome",
						draftGenomes, "GENOME");
	parser.addStoreOpt('i', "indir",
					   "directory containing input files",
					   options.indir, "DIR");
	parser.addStoreOpt('o', "outdir",
					   "directory to output result files",
					   options.outdir, "DIR");
	parser.addStoreTrueOpt(0, "output-hits",
						   "output hit files after each filtering of edges during the algorithm",
						   options.outputHits);
	parser.addStoreTrueOpt(0, "output-runs",
						   "output runs after each each clique joining",
						   options.outputRuns);
	parser.addStoreOpt(0, "pairwisehits",
					   "instead of reading from hit files, use only those hits from FILENAME",
					   options.phitFilename, "FILENAME");
	parser.addStoreTrueOpt('q', "quiet",
						   "do not output extra progress information on standard error",
						   options.quiet);
	parser.addStoreTrueOpt(0, "version",
						   "display version information and exit",
						   version);
	parser.addAppendArg("", "", nondraftGenomes);
	parser.parse(argv, argv + argc);
	
	if (version) {
		writeVersion(std::cerr);
		return EXIT_SUCCESS;
	}

	try {
		// Check that there are at least 2 genomes
		if (nondraftGenomes.size() + draftGenomes.size() < 2) {
			throw std::runtime_error("At least 2 genomes must be specified");
		}
		
		// Check that number of genomes is under limit for bitmask
		if (nondraftGenomes.size() + draftGenomes.size() > MAX_GENOMES) {
			throw std::runtime_error("At most " +
									 util::string::toString(MAX_GENOMES) +
									 " genomes may be specified");
		}
		
		// Form data directory path and check for its existence
		Path dataDir(options.indir);
		if (not dataDir.exists() or not dataDir.isDirectory()) {
			throw std::runtime_error("Invalid data directory: " + options.indir);
		}
		
		// Form output directory path and check for its existence
		Path outDir(options.outdir);
		if (not outDir.exists() or not outDir.isDirectory()) {
			throw std::runtime_error("Invalid output directory: " + options.outdir);
		}
		
		// Display parameters used
		writeVersion(std::cerr);
		std::cerr << "PARAMETERS:\n"
				  << "non-draft genomes = "
				  << util::string::join(nondraftGenomes.begin(),
										nondraftGenomes.end(), " ") << '\n'
				  << "    draft genomes = "
				  << util::string::join(draftGenomes.begin(),
										draftGenomes.end(), " ") << '\n'
				  << "            indir = " << options.indir << '\n'
				  << "           outdir = " << options.outdir << '\n'
				  << "       repeat-num = " << options.repeatNum << '\n'
				  << "       repeat-pct = " << options.repeatPct << '\n'
				  << "         max-eval = " << options.maxE << '\n'
				  << "        prune-pct = " << options.prunePct << '\n'
				  << "    join-distance = " << options.maxDist << '\n'
				  << "   min-run-length = " << options.minRunLength << '\n'
				  << "          padding = " << options.padding << '\n'
				  << '\n';
		
		// Record start time
		std::time_t startTime = std::time(0);
		
		// Initialize genomes
		for (std::vector<std::string>::iterator it = nondraftGenomes.begin();
			 it != nondraftGenomes.end(); ++it) {
			Genome::addGenome(*it, false);
		}
		// Add draft genomes at end
		for (std::vector<std::string>::iterator it = draftGenomes.begin();
			 it != draftGenomes.end(); ++it) {
			Genome::addGenome(*it, true);
		}
		
		// Load data
		cerr << "Loading input files...\n";
		Genome::loadFiles(dataDir, options.maxE, options.phitFilename);
		
		cerr << "Time spent loading files: "
			 << time(0) - startTime << " seconds\n";
	
		vector<Run*> runs;
		
		if (!options.phitFilename.empty()) {
			cerr << "Joining pairwise maps...\n";
			joinPairwiseMaps2(runs, options);
		} else {
			cerr << "Making map...\n";
			makeMap(runs, options);
		}
		
		// Output run stats
		RunStats stats = calcRunStats(runs);
		cout << "Number of runs: " << stats.numRuns << '\n'
			 << "Number of cliques: " << stats.numCliques << '\n'
			 << "Mean run length: " << stats.meanRunLength << '\n'
			 << "Median run length: " << stats.medianRunLength << '\n'
			 << "Max run length: " << stats.maxRunLength << '\n'
			 << "Min run length: " << stats.minRunLength << '\n';
		
		// Output coverage of anchors
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			Genome* genome = Genome::getGenome(g);
			size_t numCovered = genome->getNumAnchorsInRuns();
			float pctCovered = 0;
			if (genome->getNumAnchors() > 0) {
				pctCovered = 100.0 * numCovered / genome->getNumAnchors();
			}
			cout << "Coverage of " << genome->getName() << " anchors: "
				 << pctCovered << '%'
				 << " (" << numCovered << "/" << genome->getNumAnchors() << ")"
				 << '\n';
		}
		
		// Output file with genome names
		OutputFileStream genomeFile(outDir / "genomes");
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			if (g != 0) { genomeFile << '\t'; }
			genomeFile << Genome::getGenome(g)->getName();
		}
		genomeFile << '\n';
		genomeFile.close();
		
		// Output coverage files
		cerr << "Writing coverage files...\n";
		for (size_t g = 0; g < Genome::getNumGenomes(); ++g) {
			Genome* genome = Genome::getGenome(g);
			OutputFileStream covFile(outDir / (genome->getName() + ".coverage"));
			genome->printCoverage(covFile);
			cout << "Coverage of "
				 << genome->getName()
				 << ": "
				 << genome->getPctCoverage() * 100.0 << '%'
				 << '\n';
		}
		
		// Output runs
		cerr << "Writing runs...\n";
		string genomeStr = Genome::getGenomeNamesString();
		OutputFileStream runFile(outDir / "runs");
		printRuns(runs, runFile);
		
		// Output pre-map
		cerr << "Writing pre-map...\n";
		OutputFileStream preMapFile(outDir / "pre.map");
		printMap(runs, preMapFile, false);
		
		// Output extended map
		cerr << "Writing extended map...\n";
		OutputFileStream extendedMapFile(outDir / "map");
		printMap(runs, extendedMapFile, true);
		
		// Output pairwise hits
		cerr << "Writing pairwise hits...\n";
		OutputFileStream pairwiseHitFile(outDir / "pairwisehits");
		printPairwiseHits(runs, pairwiseHitFile);
		
		// Output AGP files
		cerr << "Writing AGP files...\n";
		Genome::writeAGPFiles(outDir);

		// Output anchor files
		cerr << "Writing anchor files...\n";
		Genome::writeAnchorFiles(outDir);
		
		// Output permutations
		cerr << "Writing run permutation files...\n";
		Genome::writeRunPermFiles(outDir);
		
		// Write run time
		cerr << "Run time: " << time(0) - startTime << " seconds\n";

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
