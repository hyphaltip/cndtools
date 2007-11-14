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

#include "bio/sdb.hh"
#include "util/io.hh"
#include "util/string.hh"

namespace bio { namespace SDB {

	const unsigned int DB::MAGIC_NUMBER = 0x571CD854;
	const unsigned int DB::VERSION = 0;
	
	DB::Iterator DB::begin() {
		if (!opened) {
			throw std::runtime_error("Attempted to search unopened database");
		}
		
		if (!isIndexRead()) {
			readIndex();
		} else if (!isIndexSorted()) {
			sortIndex();
		}

		return Iterator(recs.begin());
	}
		
	DB::Iterator DB::end() {
		if (!opened) {
			throw std::runtime_error("Attempted to search unopened database");
		}
		
		if (!isIndexRead()) {
			readIndex();
		} else if (!isIndexSorted()) {
			sortIndex();
		}
		
		return Iterator(recs.end());
	}
	
	DB::DB(bool cache)
		: cache(cache) {
		reset();
	}
	
	DB::~DB() {
		if (opened) {
			close();
		}
	}

	void DB::reset() {
		lastRecord = Record(this);
		opened = false;
		readonly = false;
		indexSorted = false;
		version = VERSION;
		numRecs = 0;
		seqSize = 0;
		indexSize = 0;
		recs.clear();
	}

	void DB::setCache(bool cache) {
		this->cache = cache;
	}
	
	void DB::open(const filesystem::Path& filename,
				  bool readonly,
				  bool create) {
		// Check to make sure we are not already open
		if (opened) {
			throw std::runtime_error("Attempted to open database before closing previous database");
		}

		if (readonly) {
			openForRead(filename);
		} else if (create) {
			openForNew(filename);
		} else {
			openForAppend(filename);
		}

		this->readonly = readonly;
		this->opened = true;
	}

	void DB::openForNew(const filesystem::Path& filename) {
		// Truncate any existing file
		openFile(filename, "wb");

		writeHeader();
	}

	void DB::openForAppend(const filesystem::Path& filename) {
		// Open for writing and reading so that we can read in the index
		openFile(filename, "r+b");
		
		readHeader();

		// Read entire index to memory because we will be updating it
		readIndex();

		// Position for writing where the old index began
		
		fseeko(strm, headerSize() + seqSize, SEEK_SET);
	}

	void DB::openForRead(const filesystem::Path& filename) {
		// Open for reading only
		openFile(filename, "rb");
		
		readHeader();
	}

	void DB::openFile(const filesystem::Path& filename,
					  const char* mode) {
		strm = fopen(filename.toString().c_str(), mode);
		if (!strm) {
			throw std::runtime_error(std::string("Failed to open file ") +
									 filename.toString());
		}
	}
	
	void DB::readHeader() {
		unsigned int magic;
		if (!(util::io::binary::read(strm, magic) &&
			  util::io::binary::read(strm, version) &&
			  util::io::binary::read(strm, numRecs) &&
			  util::io::binary::read(strm, seqSize) &&
			  util::io::binary::read(strm, indexSize))) {
			throw std::runtime_error("Error while reading database header");
		} else if (magic != MAGIC_NUMBER) {
			throw std::runtime_error("File is not a valid sequence database file");
		} else if (version > VERSION) {
			throw std::runtime_error("Unknown database version");
		}
	}

	void DB::writeHeader() {
		if (!(util::io::binary::write(strm, MAGIC_NUMBER) &&
			  util::io::binary::write(strm, VERSION) &&
			  util::io::binary::write(strm, numRecs) &&
			  util::io::binary::write(strm, seqSize) &&
			  util::io::binary::write(strm, indexSize))) {
			throw std::runtime_error("Error while writing database header");
		}
	}

	unsigned int DB::headerSize() const {
		return (util::io::binary::size(MAGIC_NUMBER) +
				util::io::binary::size(VERSION) +
				util::io::binary::size(numRecs) +
				util::io::binary::size(seqSize) +
				util::io::binary::size(indexSize));
	}
	
	void DB::close() {
		// Make sure we are already open
		if (!opened) {
			throw std::runtime_error("Attempted to close non-opened database");
		}
		
		// Write the index and header information if we are in write mode
		if (!readonly) {
			fseeko(strm, headerSize() + seqSize, SEEK_SET);
			writeIndex();
			fseeko(strm, 0, SEEK_SET);
			writeHeader();
		}

		fclose(strm);

		// Reset values to defaults
		reset();
	}

	unsigned int DB::calcSizes(std::vector<Record*>::iterator start,
							   std::vector<Record*>::iterator end) {
		if (start == end) {
			return 0;
		}

		std::vector<Record*>::iterator middle = start + (end - start) / 2;
		(*middle)->leftSize = calcSizes(start, middle);
		(*middle)->rightSize = calcSizes(middle + 1, end);
		return (*middle)->leftSize
			+ (*middle)->rightSize
			+ (*middle)->getBinarySize();
	}

	void DB::writeRecs(std::vector<Record*>::iterator start,
					   std::vector<Record*>::iterator end) {
		if (start == end) {
			return;
		}
		
		std::vector<Record*>::iterator middle = start + (end - start) / 2;
		(*middle)->write();
		writeRecs(start, middle);
		writeRecs(middle + 1, end);
	}

	void DB::sortIndex() {
		std::sort(recs.begin(), recs.end(), RecordSorter());
		indexSorted = true;
		for (unsigned int i = 0; i < recs.size(); ++i) {
			recs[i]->recNum = i + 1;
		}
	}
	
	void DB::writeIndex() {
		if (!isIndexSorted()) {
			sortIndex();
		}
		indexSize = calcSizes(recs.begin(), recs.end());
		writeRecs(recs.begin(), recs.end());
	}

	void DB::readIndex() {
		if (isIndexRead()) {
			return;
		}
		
		indexSorted = true;

		// Stop if there are no records
		if (numRecs == 0) {
			return;
		}
		
		// Position at start of index
		fseeko(strm, headerSize() + seqSize, SEEK_SET);
		//strm.seekg(-indexSize, std::ios::end);

		readIndexRecursive();		
	}

	void DB::readIndexRecursive() {
		Record* rec = new Record(this);
		rec->read();
		if (rec->leftSize != 0) {
			readIndexRecursive();
		}
		recs.push_back(rec);
		if (rec->rightSize != 0) {
			readIndexRecursive();
		}
	}
		
	bool DB::isIndexRead() const { return recs.size() == numRecs; }
	bool DB::isIndexSorted() const { return indexSorted; }

	std::string DB::getSeq(const std::string& title) {
		cacheRec(title);
		return lastRecord.getSeq();
	}
	
	std::string DB::getSeq(const std::string& title,
						   const unsigned int start,
						   const unsigned int end,
						   const char strand,
						   const alphabet::Nucleotide& alphabet) {
		cacheRec(title);
		return lastRecord.getSeq(start, end, strand, alphabet);
	}

	std::string DB::getSeq(const genome::Interval& i,
					   const alphabet::Nucleotide& alphabet) {
		return getSeq(i.getChrom(),
					  i.getStart(),
					  i.getEnd(),
					  i.getStrand(),
					  alphabet);
	}
	
	std::string DB::getSeq(const unsigned int recNum) {
		cacheRec(recNum);
		return lastRecord.getSeq();
	}
	
	std::string DB::getSeq(const unsigned int recNum,
						   const unsigned int start,
						   const unsigned int end,
						   const char strand,
						   const alphabet::Nucleotide& alphabet) {
		cacheRec(recNum);
		return lastRecord.getSeq(start, end, strand, alphabet);
	}

	void DB::getRec(const std::string& title,
					Record& rec) {
		cacheRec(title);
		// Avoid copying and disturbing caching if rec is already correct
		if (!(rec == lastRecord)) {
			rec = lastRecord;
		}
	}

	void DB::getRec(const unsigned int recNum,
					Record& rec) {
		cacheRec(recNum);
		// Avoid copying and disturbing caching if rec is already correct
		if (!(rec == lastRecord)) {
			rec = lastRecord;
		}
	}

	void DB::cacheRec(const std::string& title) {
		// Check to make sure database has been opened
		if (!opened) {
			throw std::runtime_error("Attempted to search unopened database");
		}
		// Check if this record is already stored
		if (!lastRecord.getTitle().empty() && lastRecord.getTitle() == title) {
			return;
		}
		// Check if there are no records in the database
		if (numRecs == 0 || !lookupRec(title)) {
			throw NotFoundError("Record \"" + title + "\" not in database");
		}
	}
	
	void DB::cacheRec(const unsigned int recNum) {
		// Check to make sure database has been opened
		if (!opened) {
			throw std::runtime_error("Attempted to search unopened database");
		}
		// Check if this record is already stored
		if (!lastRecord.getTitle().empty() && lastRecord.getRecNum() == recNum) {
			return;
		}
		// Check if there are no records in the database
		if (numRecs == 0 || !lookupRec(recNum)) {
			throw NotFoundError("Record \"" + util::string::toString(recNum)
								+ "\" not in database");
		}
	}

	bool DB::lookupRec(const std::string& title) {
		// Use in-memory index if possible
		if (isIndexRead()) {
			if (!isIndexSorted()) {
				sortIndex();
			}
			return lookupRec(title, recs.begin(), recs.end());
		}
		// Else, look record up in disk index
		else {			
			// Position at the start of the index
			fseeko(strm, headerSize() + seqSize, SEEK_SET);
			
			while (true) {
				lastRecord.read();
				int cmp = title.compare(lastRecord.title);
				if (cmp == 0) {
					return true;
				} else if (cmp > 0) {
					if (lastRecord.rightSize == 0) {
						return false;
					} else {
						// Skip over the left subtree and to the start of
						// the right subtree
						fseeko(strm, lastRecord.leftSize, SEEK_CUR);
						continue;
					}
				} else {
					if (lastRecord.leftSize == 0) {
						return false;
					} else {
						// The read position should already be at the
						// beginning of the left subtree
						continue;
					}
				}
			}
		}
	}	

	bool DB::lookupRec(const unsigned int recNum) {
		// Use in-memory index if possible
		if (isIndexRead()) {
			if (!isIndexSorted()) {
				sortIndex();
			}
			return lookupRec(recNum, recs.begin(), recs.end());
		}
		// Else, look record up in disk index
		else {
			// Position at the start of the index
			fseeko(strm, headerSize() + seqSize, SEEK_SET);
			
			while (true) {
				lastRecord.read();
				if (recNum == lastRecord.getRecNum()) {
					return true;
				} else if (recNum > lastRecord.getRecNum()) {
					if (lastRecord.rightSize == 0) {
						return false;
					} else {
						// Skip over the left subtree and to the start of
						// the right subtree
						fseeko(strm, lastRecord.leftSize, SEEK_CUR);
						continue;
					}
				} else {
					if (lastRecord.leftSize == 0) {
						return false;
					} else {
						// The read position should already be at the
						// beginning of the left subtree
						continue;
					}
				}
			}
		}
	}	
	
	bool DB::lookupRec(const std::string& title,
					   std::vector<Record*>::iterator begin,
					   std::vector<Record*>::iterator end) {
		if (begin == end) {
			return false;
		}
		std::vector<Record*>::iterator middle = begin + (end - begin) / 2;
		int cmp = title.compare((*middle)->title);
		if (cmp == 0) {
			lastRecord = **middle;
			return true;
		} else if (cmp > 0) {
			return lookupRec(title, middle + 1, end);
		} else {
			return lookupRec(title, begin, middle);
		}
	}

	bool DB::lookupRec(const unsigned int recNum,
					   std::vector<Record*>::iterator begin,
					   std::vector<Record*>::iterator end) {
		if (begin == end) {
			return false;
		}
		std::vector<Record*>::iterator middle = begin + (end - begin) / 2;
		if (recNum == (*middle)->recNum) {
			lastRecord = **middle;
			return true;
		} else if (recNum > (*middle)->recNum) {
			return lookupRec(recNum, middle + 1, end);
		} else {
			return lookupRec(recNum, begin, middle);
		}
	}
	
	void DB::putRec(const std::string& title,
					const std::string& sequence,
					const bool compressed) {
		if (!opened) {
			throw std::runtime_error("Attempted to write to unopened database");
		} else if (readonly) {
			throw std::runtime_error("Attempted to write database opened for reading");
		} else if (title.empty()) {
			throw std::runtime_error("Attempted to write record with empty title");
		}

		// Seek to the end of the stored sequences
		fseeko(strm, headerSize() + seqSize, SEEK_SET);

		Record* rec = new Record(this);
		rec->title = title;
		rec->seqLength = sequence.length();
		rec->seqPos = headerSize() + seqSize;
		rec->compressed = (compressed ? 1 : 0);
		recs.push_back(rec);
		++numRecs;
		if (compressed) {
			std::string encoded = alphabet::Nucleotide::nibEncode(sequence);
			util::io::binary::write(strm, encoded.data(), encoded.length());
			seqSize += encoded.length();
		} else {
			util::io::binary::write(strm, sequence.data(), sequence.length());
			seqSize += sequence.length();
		}
		indexSorted = false;
	}

	// Implementation of Record

	Record::Record()
		: db(NULL) {}
	
	Record::Record(DB* db)
		: db(db) {}

	bool Record::operator==(const Record& other) const {
		return db == other.db
			&& recNum == other.recNum
			&& seqPos == other.seqPos
			&& seqLength == other.seqLength
			&& compressed == other.compressed;
	}
	
	void Record::write() const {
		if (!(util::io::binary::write(db->strm, title)
			  && util::io::binary::write(db->strm, recNum)
			  && util::io::binary::write(db->strm, leftSize)
			  && util::io::binary::write(db->strm, rightSize)
			  && util::io::binary::write(db->strm, seqLength)
			  && util::io::binary::write(db->strm, seqPos)
			  && util::io::binary::write(db->strm, compressed))) {
			throw std::runtime_error("Error while reading sdb record");
		}
	}

	void Record::read() {
		if (!(util::io::binary::read(db->strm, title)
			  && util::io::binary::read(db->strm, recNum)
			  && util::io::binary::read(db->strm, leftSize)
			  && util::io::binary::read(db->strm, rightSize)
			  && util::io::binary::read(db->strm, seqLength)
			  && util::io::binary::read(db->strm, seqPos)
			  && util::io::binary::read(db->strm, compressed))) {
			throw std::runtime_error("Error while reading sdb record");
		}
		// Clear sequence cache
		seq.clear();
	}
	
	unsigned int Record::getBinarySize() const {
		return util::io::binary::size(title)
			+ util::io::binary::size(recNum)
			+ util::io::binary::size(leftSize)
			+ util::io::binary::size(rightSize)
			+ util::io::binary::size(seqLength)
			+ util::io::binary::size(seqPos)
			+ util::io::binary::size(compressed);
	}

	const std::string& Record::getTitle() const { return title; }
	unsigned int Record::getRecNum() const { return recNum; }
	unsigned int Record::getLength() const { return seqLength; }

	std::string Record::getSeq() {
		return getSeq(0, seqLength);
	}

	std::string Record::getSeq(const unsigned int start,
							   const unsigned int end,
							   const char strand,
							   const alphabet::Nucleotide& alphabet) {
		if (start < 0 || end > seqLength || end < start) {
			throw OutOfBoundsError("Bad coordinates for " + title
								   + ": " + util::string::toString(start)
								   + "-" + util::string::toString(end));
		}
		
		std::string seq;

		if (start == end) {
			return seq;
		} else if (compressed) {
			unsigned int startByte = start / 2;
			unsigned int offset = start % 2;
			unsigned int endByte = (start == end ? startByte : (end + 1) / 2);
			
			std::string encoded = readSeq(startByte, endByte - startByte);
			seq = alphabet::Nucleotide::nibDecode(encoded).substr(offset, end - start);
		} else {
			seq = readSeq(start, end - start);
		}
		
		if (strand == '-') {
			alphabet.reverseComplementInPlace(seq);
		}
		
		return seq;
	}

	size_t Record::getCompressedLength() const {
		if (compressed) {
			return (seqLength / 2) + (seqLength % 2);
		} else {
			return seqLength;
		}
	}

	std::string Record::readSeq(const unsigned int offset,
								const unsigned int length) {
		if (db->cache) {
			if (seq.empty() && seqLength != 0) {
				std::vector<char> buffer(getCompressedLength());
				bufferSeq(0, buffer);
				seq.assign(buffer.begin(), buffer.end());
			}

			if (offset == 0 && length == seq.length()) {
				return seq;
			} else {
				return seq.substr(offset, length);
			}
			
		} else {
			std::vector<char> buffer(length);
			bufferSeq(offset, buffer);
			return std::string(buffer.begin(), buffer.end());
		}
	}

	void Record::bufferSeq(const unsigned int offset,
						   std::vector<char>& buffer) {
		fseeko(db->strm, seqPos + offset, SEEK_SET);
		if (!util::io::binary::read(db->strm, &buffer[0], buffer.size())) {
			throw std::runtime_error("Error while reading sdb sequence");
		}
	}
	
} }
