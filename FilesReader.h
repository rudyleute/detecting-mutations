#ifndef FILESREADER_H
#define FILESREADER_H
#include <string>
#include <set>
#include <boost/icl/interval_map.hpp>
#include <utility>

using namespace std;
using namespace boost::icl;

struct NucleoCounter;
struct AlignmentMaps;


class FilesReader {
public:
	static std::map<char, size_t> nucleoMapping;
	static std::map<size_t, char> rNucleoMapping;
	static std::map<string, size_t> cigarIndices;

	static AlignmentMaps getAlignments(const std::string& fileName, const size_t &from, const size_t &to);
	static size_t getReferenceLength(const std::string& fileName);
	static char complementSymbol(const char &n);
	static void reverseComplement(string &read);
	static string getRefGen(const string &fileName);
	static std::map<size_t, std::pair<char, char>> readReferenceCsv(const string &fileName);
};

struct NucleoCounter {
private:
	std::vector<size_t> counters;

public:
	NucleoCounter(): counters(FilesReader::nucleoMapping.size()) {}

	void increase(const char &nucleo) {
		counters[FilesReader::nucleoMapping[nucleo]]++;
	}

	void setCounter(const char &nucleo, const size_t &value) {
		counters[FilesReader::nucleoMapping[nucleo]] = value;
	}

	char findMax(const char &base, const bool isInsertion=false) const {
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
		auto result = curMaxNucleos.find(base);

		//If there are more than two nucleotides with the same score and one of them is the base
		if (isInsertion && curMaxNucleos.size() > 1 && result != curMaxNucleos.end()) {
			for (auto nucleo: curMaxNucleos) {
				if (nucleo != base) return nucleo;
			}
		};
		if (curMaxNucleos.find(base) == curMaxNucleos.end()) return *curMaxNucleos.begin();
		return base;
	}

	void flush() {
		std::fill(counters.begin(), counters.end(), 0);
	}

	size_t size() const {
		size_t counter = 0;
		for (const size_t i : counters) counter += i;

		return counter;
	}
};

struct AlignmentMaps {
	std::map<size_t, set<pair<string, string>>> startingPos;
	std::map<size_t, pair<NucleoCounter, set<string>>> insertionCount;

	AlignmentMaps(std::map<size_t, set<pair<string, string>>> alignments, std::map<size_t, pair<NucleoCounter, set<string>>> insertions): startingPos(move(alignments)), insertionCount(move(insertions)) {}
};

#endif //FILESREADER_H
