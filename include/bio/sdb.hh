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

#ifndef __BIO_SDB_HH__
#define __BIO_SDB_HH__

#include <stdexcept>
#include <fstream>

#include "boost/iterator/indirect_iterator.hpp"

#include "bio/genome/Interval.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "bio/alphabet/AmbiguousNucleotide.hh"
#include "filesystem/Path.hh"

namespace bio {

namespace SDB {

	class NotFoundError : public std::runtime_error {
	public:
		NotFoundError(const std::string& what);
	};

	class OutOfBoundsError : public std::runtime_error {
	public:
		OutOfBoundsError(const std::string& what);
	};
	
	inline NotFoundError::NotFoundError(const std::string& what)
		: std::runtime_error(what) {
	}

	inline OutOfBoundsError::OutOfBoundsError(const std::string& what)
		: std::runtime_error(what) {
	}
	
	class DB;
	
	class Record {
	public:
		Record();
		
		const std::string& getTitle() const;
		unsigned int getRecNum() const;
		unsigned int getLength() const;
		std::string getSeq();
		std::string getSeq(const unsigned int start,
						   const unsigned int end,
						   const char strand='+',
						   const alphabet::Nucleotide& alphabet=alphabet::AmbiguousDNA);

		bool operator==(const Record& other) const;
		
	private:
		Record(DB* db);
		void write() const;
		void read();
		unsigned int getBinarySize() const;
		std::string readSeq(const unsigned int offset,
							const unsigned int length);
		void bufferSeq(const unsigned int offset,
					   std::vector<char>& buffer);
		size_t getCompressedLength() const;

		DB* db;
		std::string seq;
		
		std::string title;
		unsigned int recNum;
		unsigned int leftSize;
		unsigned int rightSize;
		unsigned int seqLength;
		off_t seqPos;
		unsigned char compressed;

		friend class DB;
		friend class RecordSorter;
	};

	struct RecordSorter {
		bool operator()(const Record* r1, const Record* r2) const {
			return r1->title < r2->title;
		}
	};
	
	class DB {
	public:
		typedef boost::indirect_iterator<std::vector<Record*>::iterator> Iterator;

		DB(bool cache=false);
		~DB();

		void open(const filesystem::Path& filename,
				  bool readonly=true,
				  bool create=false);

		void close();

		void putRec(const std::string& title,
					const std::string& sequence,
					const bool compressed=false);
		
		void getRec(const std::string& title, Record& rec);
		void getRec(const unsigned int recNum, Record& rec);

		std::string getSeq(const std::string& title);
		std::string getSeq(const std::string& title,
						   const unsigned int start,
						   const unsigned int end,
						   const char strand='+',
						   const alphabet::Nucleotide& alphabet=alphabet::AmbiguousDNA);

		std::string getSeq(const genome::Interval& i,
						   const alphabet::Nucleotide& alphabet=alphabet::AmbiguousDNA);

		std::string getSeq(const unsigned int recNum);
		std::string getSeq(const unsigned int recNum,
						   const unsigned int start,
						   const unsigned int end,
						   const char strand='+',
						   const alphabet::Nucleotide& alphabet=alphabet::AmbiguousDNA);
		
		Iterator begin();
		Iterator end();

		void setCache(bool cache=true);
		void readIndex();

	private:
		void openFile(const filesystem::Path& filename,
					  const char* openmode);
		void openForNew(const filesystem::Path& filename);
		void openForAppend(const filesystem::Path& filename);
		void openForRead(const filesystem::Path& filename);

		void readHeader();
		void writeHeader();
		unsigned int headerSize() const;

		unsigned int calcSizes(std::vector<Record*>::iterator start,
							   std::vector<Record*>::iterator end);
		void writeRecs(std::vector<Record*>::iterator start,
					   std::vector<Record*>::iterator end);

		void writeIndex();
		void readIndexRecursive();
		void sortIndex();

		bool isIndexRead() const;
		bool isIndexSorted() const;

		void cacheRec(const std::string& title);
		void cacheRec(const unsigned int recNum);
		
		bool lookupRec(const std::string& title);
		bool lookupRec(const unsigned int recNum);
		
		bool lookupRec(const std::string& title,
					   std::vector<Record*>::iterator begin,
					   std::vector<Record*>::iterator end);
		bool lookupRec(const unsigned int recNum,
					   std::vector<Record*>::iterator begin,
					   std::vector<Record*>::iterator end);

		bool lookupRecInBuffer(const std::string& title,
							   std::istream& buffer);
		
		void reset();

		bool cache;
		Record lastRecord;
		
		bool opened;
		bool readonly;
		bool indexSorted;
		FILE* strm;
		unsigned int version;
		unsigned int numRecs;
		off_t seqSize;
		unsigned int indexSize;

		std::vector<Record*> recs;
	
		static const unsigned int MAGIC_NUMBER;
		static const unsigned int VERSION;

		friend class Record;
	};
	
};

};

#endif // __BIO_SDB_HH__
