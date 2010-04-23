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

#ifndef __GENOME_HH__
#define __GENOME_HH__

#include <string>

#include "types.hh"
#include "bio/sdb.hh"
#include "boost/unordered_map.hpp"
#include "filesystem/Path.hh"

typedef boost::unordered_map<std::string, Genome*> GenomeMap;
typedef std::vector<Genome*> GenomeList;

class Genome {
private:
	static filesystem::Path dbDir;
	static size_t numGenomes;
	static GenomeMap genomeMap;
	static GenomeList genomeList;
	
	std::string name;
	size_t num;
	bio::SDB::DB db;
public:
	Genome(const std::string& name);

	std::string getName() const;
	size_t getNum() const;
	size_t getChromLen(const std::string& chrom);
			
	std::string getSeq(const std::string& chrom,
					   const size_t start,
					   const size_t end,
					   const char strand='+');
	
	static void setDBDir(const std::string& dir);
	static size_t getNumGenomes();
	static Genome* getGenome(const std::string& name, bool create = false);
	static Genome* getGenome(size_t num);
};

#endif // __GENOME_HH__
