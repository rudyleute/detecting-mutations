#ifndef FILESREADER_H
#define FILESREADER_H
#include <string>
#include <set>
#include <boost/icl/interval_map.hpp>
#include <utility>
#include <list>
#include <htslib/sam.h>
#include <iostream>

using namespace std;
using namespace boost::icl;


struct NucleoCounter;
struct AlignmentMaps;
struct Insertions;

using InsertionMap = std::map<size_t, std::pair<NucleoCounter, std::set<std::string>>>;
using CigarString = std::list<pair<char, size_t>>;
using Mutations = std::map<size_t, std::pair<char, char>>;


class FilesReader {
	static Mutations getVCFInsertions(const string& ref, const string& alt, const size_t& pos);

public:
	static std::map<char, size_t> nucleoMapping;
	static std::map<size_t, char> rNucleoMapping;
	static std::map<std::string, std::tuple<size_t, size_t, size_t>> cigarIndices;

	static AlignmentMaps getAlignments(
		const std::string& fileName,
		const size_t& from,
		const size_t& to,
		InsertionMap prevIterInsertions
	);
	static size_t getReferenceLength(const std::string& fileName);
	static string getRefGen(const string& fileName);
	static std::map<size_t, std::pair<char, char>> readFreeBayesVCF(const string& fileName);
	static CigarString getCigarString(const bam1_t* b);
	static string getRead(const bam1_t* b);
	static string getExpandedRead(string read, CigarString& cigar);
	static string formFullPath(const string& fileName);
};

struct NucleoCounter {
private:
	std::vector<size_t> counters;

	vector<size_t> getCounters() const {
		return counters;
	}

public:
	NucleoCounter(): counters(FilesReader::nucleoMapping.size()) {}

	void increase(const char& nucleo) {
		counters[FilesReader::nucleoMapping[nucleo]]++;
	}

	void setCounter(const char& nucleo, const size_t& value) {
		counters[FilesReader::nucleoMapping[nucleo]] = value;
	}

	char findMax(const char& base) const {
		size_t curMax = 0;
		set<char> curMaxNucleos;

		for (size_t i = 0; i != counters.size(); i++) {
			if (counters[i] >= curMax) {
				if (counters[i] > curMax) {
					curMax = counters[i];
					curMaxNucleos.clear();
				}

				curMaxNucleos.insert(FilesReader::rNucleoMapping[i]);
			}
		}

		const double ratio = static_cast<double>(curMax) / size();
		//TODO probably it is worth considering adding an epsilon here?
		if (ratio >= 0.5) {
			//Return the first alphabetically sorted non-base value if it is seen in at least 50% of cases
			if (ratio > 0.5 || curMaxNucleos.size() == 1) return *curMaxNucleos.begin();
			for (const auto& elem : curMaxNucleos) if (elem != base) return elem;
		}

		return base;
	}

	void flush() {
		std::fill(counters.begin(), counters.end(), 0);
	}

	void merge(const NucleoCounter& other) {
		const auto result = other.getCounters();
		for (size_t i = 0; i != counters.size(); i++) {
			counters[i] += result[i];
		}
	}

	size_t size() const {
		size_t counter = 0;
		for (const size_t i : counters) counter += i;

		return counter;
	}
};

struct Insertions {
private:
	InsertionMap insertions;
	InsertionMap nextWindowInsertions;
	set<size_t> insertionIndices;

	InsertionMap* curMap;
	string expandedRead;
	string name;

	void addValues(const size_t& refGenIndex, const size_t& curReadIndex, const size_t& start, const size_t& end, const bool isInsertion) {
		for (size_t i = start; i != end; i++) {
			if (isInsertion) {
				(*curMap)[refGenIndex + i].first.increase(expandedRead[curReadIndex + i]);
				(*curMap)[refGenIndex + i].second.insert(name);
				insertionIndices.insert(refGenIndex + i);
			} else if ((*curMap)[refGenIndex + i].second.find(name) == (*curMap)[refGenIndex + i].second.end()) {
				(*curMap)[refGenIndex + i].first.increase('-');
			}
		}
	}

public:
	Insertions() = default;

	Insertions(
		const InsertionMap& insertions,
		const InsertionMap& nextWindowInsertions = InsertionMap(),
		const set<size_t>& insertionIndices = set<size_t>()
	): insertions(insertions), nextWindowInsertions(nextWindowInsertions), insertionIndices(insertionIndices) {
		curMap = &this->insertions;
	}

	void setRead(const string& expandedRead, const string& name) {
		this->expandedRead = move(expandedRead);
		this->name = move(name);
	}
	InsertionMap getNextWindowInsertions() const {
		return nextWindowInsertions;
	}

	void addInsertion(
		const size_t& refGenIndex,
		const size_t& curReadIndex,
		const size_t& end,
		const size_t& left,
		const bool isInsertion = false
	) {
		if (left != end) {
			this->addValues(refGenIndex, curReadIndex, 0, left, isInsertion);
			curMap = &this->nextWindowInsertions;
			this->addValues(refGenIndex, curReadIndex, left, end, isInsertion);
			curMap = &this->insertions;
		} else this->addValues(refGenIndex, curReadIndex, 0, end, isInsertion);
	}

	Mutations findInsertionMutations(const size_t &minReads) {
		Mutations errors;

		for (const size_t &index: insertionIndices) {
			if (insertions[index].first.size() >= minReads)
				if (const char maxNucleo = insertions[index].first.findMax('-'); maxNucleo != '-') errors[index] = make_pair(
					maxNucleo, 'I');
		}

		return errors;
	}
};

struct AlignmentMaps {
	std::map<size_t, set<pair<string, string>>> startingPos;
	Insertions windowInsertions;

	AlignmentMaps(
		const std::map<size_t, set<pair<string, string>>>& alignments,
		const Insertions& windowInsertions
	): startingPos(alignments),
	   windowInsertions(windowInsertions) {}

	AlignmentMaps() = default;
};

#endif //FILESREADER_H
