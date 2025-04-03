#ifndef SAMREADER_H
#define SAMREADER_H
#include <string>
#include <boost/icl/interval_map.hpp>
#include <utility>

using namespace std;
using namespace boost::icl;

struct AlignmentMaps {
	interval_map<size_t, std::string> intervalAlignment;
	std::map<size_t, size_t> insertionCount;

	AlignmentMaps(interval_map<size_t, std::string> alignments, std::map<size_t, size_t> insertions): intervalAlignment(std::move(alignments)), insertionCount(std::move(insertions)) {}
};

class SamReader {
public:
	static AlignmentMaps getAlignments(const std::string& fileName);
	static char complementSymbol(const char &n);
	static void reverseComplement(std::string &read);
};

#endif //SAMREADER_H
