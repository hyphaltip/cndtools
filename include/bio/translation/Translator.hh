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

#ifndef __BIO_TRANSLATION_TRANSLATOR_HH__
#define __BIO_TRANSLATION_TRANSLATOR_HH__

#include "bio/translation/Table.hh"

namespace bio { namespace translation {

    class Translator {
	public:
		Translator(size_t tableNum);
		Translator(const Table* table);

		std::string operator()(const std::string& seq) const;
		
		std::string translate(const std::string& seq,
							  const unsigned int phase=0,
							  const bool toStop=false) const;
		
		char translateCodon(const char base1,
							const char base2,
							const char base3) const;

		char translateCodon(const std::string& codon) const;
		
		bool isStartCodon(const char base1,
						  const char base2,
						  const char base3) const;

		bool isStartCodon(const std::string& codon) const;

		bool isStopCodon(const char base1,
						 const char base2,
						 const char base3) const;

		bool isStopCodon(const std::string& codon) const;

	private:
		const Table* table;
	};

} }

#endif // __BIO_TRANSLATION_TRANSLATOR_HH__
