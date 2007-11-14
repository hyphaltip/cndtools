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

#include "Genome.hh"

Genome::Genome(const std::string& name)
	: name(name),
	  num(numGenomes++),
	  db()
{
	filesystem::Path dbFilepath = dbDir / (name + ".sdb");
	db.open(dbFilepath);
	genomeMap.insert(std::make_pair(name, this));
	genomeList.push_back(this);
}

std::string Genome::getName() const {
	return name;
}

size_t Genome::getNum() const { return num; }

size_t Genome::getChromLen(const std::string& chrom) {
	bio::SDB::Record rec;
	db.getRec(chrom, rec);
	return rec.getLength();
}
			
std::string Genome::getSeq(const std::string& chrom,
						   const size_t start,
						   const size_t end,
						   const char strand) {
	return db.getSeq(chrom, start, end, strand);
}

filesystem::Path Genome::dbDir = std::string(".");
size_t Genome::numGenomes = 0;
GenomeMap Genome::genomeMap = GenomeMap();
GenomeList Genome::genomeList = GenomeList();

void Genome::setDBDir(const std::string& dir) { dbDir = dir; }

size_t Genome::getNumGenomes() { return numGenomes; }

Genome* Genome::getGenome(const std::string& name, bool create) {
	GenomeMap::iterator it = genomeMap.find(name);
	if (it == genomeMap.end()) {
		return create ? new Genome(name) : NULL;
	} else {
		return it->second;
	}
}
Genome* Genome::getGenome(size_t num) { return genomeList[num]; }
