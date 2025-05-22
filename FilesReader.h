#ifndef FILESREADER_H
#define FILESREADER_H
#include <string>
#include <set>
#include <boost/icl/interval_map.hpp>
#include <utility>
#include <list>
#include <htslib/sam.h>

using namespace std;
using namespace boost::icl;


struct NucleoCounter;
struct AlignmentMaps;
using InsertionMap = std::map<size_t, std::pair<NucleoCounter, std::set<std::string>>>;
using CigarString = std::list<pair<char, size_t>>;
using Mutations = std::map<size_t, std::pair<char, char>>;


class FilesReader {
	static Mutations getVCFInsertions(const string& ref, const string& alt, const size_t& pos);

public:
	static std::map<char, size_t> nucleoMapping;
	static std::map<size_t, char> rNucleoMapping;
	static std::map<std::string, std::tuple<size_t, size_t, size_t>> cigarIndices;

	static AlignmentMaps getAlignments(const std::string& fileName, const size_t& from, const size_t& to, InsertionMap prevIterInsertions);
	static size_t getReferenceLength(const std::string& fileName);
	static string getRefGen(const string& fileName);
	static std::map<size_t, std::pair<char, char>> readFreeBayesVCF(const string& fileName);
	static CigarString getCigarString(const bam1_t* b);
	static string getRead(const bam1_t* b);
	static string getExpandedRead(string read, CigarString &cigar);
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
			for (const auto &elem: curMaxNucleos)
				if (elem != base) return elem;
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

struct AlignmentMaps {
	std::map<size_t, set<pair<string, string>>> startingPos;
	InsertionMap insertions;
	InsertionMap insertionsOOB;

	AlignmentMaps(const std::map<size_t, set<pair<string, string>>>& alignments,
	              const InsertionMap& insertions,
	              const InsertionMap& insertionsOOB): startingPos(alignments),
	                                                  insertions(insertions),
	                                                  insertionsOOB(insertionsOOB) {}

	AlignmentMaps() = default;
};

#endif //FILESREADER_H
