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

#include <string>
#include <vector>

#include "util/charmanip.hh"

#include "bio/translation/Table.hh"

namespace bio { namespace translation {

	Table::Table(const char* table,
				 unsigned int num,
				 std::string name,
				 std::string name2) :
		table(table),
		num(num),
		name(name),
		name2(name2)
	{
		if (num >= tables.size()) {
			tables.resize(num + 1, NULL);
		}
		tables[num] = this;
	}
	
	const Table* Table::getTable(unsigned int num) {
		return num >= tables.size() ? NULL : tables[num];
	}

	unsigned int Table::getNum() const { return num; }
	std::string Table::getName() const { return name; }
	std::string Table::getName2() const { return name2; }
		
	char Table::translate(const char base1,
						  const char base2,
						  const char base3) const {
		// 			bool masked = std::islower(base1) ||
		// 				std::islower(base2) || std::islower(base3);
			
		char aa = table[(ntEncoder(base1) << 8) +
						(ntEncoder(base2) << 4) +
						ntEncoder(base3)];
		//			return masked ? std::tolower(aa) : aa;
		return aa;
	}

	std::string Table::translate(const std::string& seq,
								 const unsigned int phase,
								 const bool toStop) const {
		unsigned int proteinSize = (seq.size() < phase ?
									0 : (seq.size() - phase) / 3);
		std::string protein(proteinSize, '?');
		for (unsigned int i = 0; i < proteinSize; ++i) {
			unsigned int start = i * 3 + phase;
			char aa = translate(seq[start],
								seq[start + 1],
								seq[start + 2]);
			if (toStop && aa == '*') {
				protein.resize(i);
				break;
			} else {
				protein[i] = aa;
			}
		}
		return protein;
	}

	bool Table::isStartCodon(const std::string& codon) const {
		return isStartCodon(codon[0], codon[1], codon[2]);
	}

	bool Table::isStopCodon(const std::string& codon) const {
		return isStopCodon(codon[0], codon[1], codon[2]);
	}
	
	bool Table::isStartCodon(const char base1,
							 const char base2,
							 const char base3) const {
		char first = std::toupper(base1);
		if (num == 1) {
			return std::toupper(base2) == 'T'
				&& std::toupper(base3) == 'G'
				&& (first == 'A' || first == 'C' || first == 'T'
					|| first == 'Y' || first == 'M' || first == 'W'
					|| first == 'H');
		} else {
			abort();
		}
		return false;
	}

	bool Table::isStopCodon(const char base1,
							const char base2,
							const char base3) const {
		return translate(base1, base2, base3) == '*';
	}
		

	const util::charmanip::Encoder Table::ntEncoder("ACGTUNMRWSYKVHDB", false, 'N');
	std::vector<const Table*> Table::tables = std::vector<const Table*>();

	#include "translation.tables"

} }
