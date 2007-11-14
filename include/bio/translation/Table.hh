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

#ifndef __BIO_TRANSLATION_TABLE_HH__
#define __BIO_TRANSLATION_TABLE_HH__

#include <string>
#include <vector>

#include "util/charmanip.hh"

namespace bio { namespace translation {

    class Table {
	public:
		Table(const char* table,
			  unsigned int num,
			  std::string name,
			  std::string name2="");
		
		static const Table* getTable(unsigned int num);

		unsigned int getNum() const;
		std::string getName() const;
		std::string getName2() const;
		
		char translate(const char base1,
					   const char base2,
					   const char base3) const;

		std::string translate(const std::string& seq,
							  const unsigned int phase=0,
							  const bool toStop=false) const;

		bool isStartCodon(const std::string& codon) const;
		bool isStopCodon(const std::string& codon) const;
		
		bool isStartCodon(const char base1,
						  const char base2,
						  const char base3) const;

		bool isStopCodon(const char base1,
						 const char base2,
						 const char base3) const;
	private:
		static const util::charmanip::Encoder ntEncoder;
		static std::vector<const Table*> tables;
		
		const char* table;
		unsigned int num;
		std::string name;
		std::string name2;
	};
	
} }

#endif // __BIO_TRANSLATION_TABLE_HH__
