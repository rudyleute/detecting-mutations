#ifndef FILESREADER_H
#define FILESREADER_H
#include <string>
#include <boost/icl/interval_map.hpp>
#include <htslib/sam.h>

#include "Structures.h"

using namespace std;
using namespace boost::icl;

class FilesManipulator {
	static MutationsVCF getVCFInsertions(const string& ref, const string& alt, const size_t& pos);

public:
	static std::map<std::string, std::tuple<size_t, size_t, size_t>> cigarIndices;

	static AlignmentMaps getAlignments(
		const string& fileName,
		const size_t& from,
		const size_t& to,
		const string& refName,
		const InsertionMap& prevIterInsertions
	);
	static size_t getRefGenLength(const std::string& fileName);
	static string getRefGenName(const string& fileName);
	static string getRefGen(const string& fileName);
	static MutationsVCF readFreeBayesVCF(const string& fileName);
	static CigarString getCigarString(const bam1_t* b);
	static string getRead(const bam1_t* b);
	static string getExpandedRead(string read, CigarString& cigar);
	static string formFullPath(const string& fileName);
	static void saveToCsv(const string& geneName, CompRes& errors);
};
#endif //FILESREADER_H
