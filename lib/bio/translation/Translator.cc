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

#include "bio/translation/Translator.hh"
#include "util/string.hh"
using util::string::toString;

namespace bio { namespace translation {

	Translator::Translator(size_t tableNum)
		: table(Table::getTable(tableNum)) {
		if (table == NULL) {
			throw std::runtime_error("No translation table numbered " +
									 toString(tableNum));
		}
	}
	
	Translator::Translator(const Table* table) : table(table) {}
	
	std::string Translator::operator()(const std::string& seq) const {
		return translate(seq);
	}
	
	std::string Translator::translate(const std::string& seq,
									  const unsigned int phase,
									  const bool toStop) const {
		return table->translate(seq, phase, toStop);
	}
	
	char Translator::translateCodon(const char base1,
									const char base2,
									const char base3) const {
		return table->translate(base1, base2, base3);
	}
	
	char Translator::translateCodon(const std::string& codon) const {
		assert(codon.size() == 3);
		return table->translate(codon[0], codon[1], codon[2]);
	}
	
	bool Translator::isStartCodon(const char base1,
								  const char base2,
								  const char base3) const {
		return table->isStartCodon(base1, base2, base3);
	}
	
	bool Translator::isStartCodon(const std::string& codon) const {
		return table->isStartCodon(codon);
	}
	
	bool Translator::isStopCodon(const char base1,
								 const char base2,
								 const char base3) const {
		return table->isStopCodon(base1, base2, base3);
	}
	
	bool Translator::isStopCodon(const std::string& codon) const {
		return table->isStopCodon(codon);
	}
	
} }
