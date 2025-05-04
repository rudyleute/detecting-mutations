#ifndef FILESREADER_H
#define FILESREADER_H
#include <string>
#include <set>
#include <boost/icl/interval_map.hpp>
#include <utility>

using namespace std;
using namespace boost::icl;

struct AlignmentMaps {
	std::map<size_t, set<string>> startingPos;
	std::map<size_t, std::map<char, int>> insertionCount;

	AlignmentMaps(std::map<size_t, set<string>> alignments, std::map<size_t, std::map<char, int>> insertions): startingPos(move(alignments)), insertionCount(move(insertions)) {}
};

class FilesReader {
public:
	static AlignmentMaps getAlignments(const string& fileName);
	static char complementSymbol(const char &n);
	static void reverseComplement(string &read);
	static string getRefGen(const string &fileName);
	static std::map<size_t, std::pair<char, char>> readReferenceCsv(const string &fileName);
};

#endif //FILESREADER_H
