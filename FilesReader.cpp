#include <iostream>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include "FilesReader.h"
#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval.hpp>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <string>

using namespace BamTools;
using namespace std;

char FilesReader::complementSymbol(const char& n) {
	switch (n) {
		case 'A': return 'T';
		case 'T': return 'A';
		case 'G': return 'C';
		case 'C': return 'G';
		case '-': return '-';
		default: return 'N';
	}
}

std::map<char, size_t> FilesReader::nucleoMapping = {
	{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}, {'-', 4}
};

std::map<size_t, char> FilesReader::rNucleoMapping = [] {
	std::map<size_t, char> reverse;
	for (const auto& [key, value] : FilesReader::nucleoMapping) {
		reverse[value] = key;
	}
	return reverse;
}();

void FilesReader::reverseComplement(std::string& read) {
	reverse(read.begin(), read.end());
	transform(read.begin(), read.end(), read.begin(), complementSymbol);
}

string FilesReader::getRefGen(const string& fileName) {
	string refGen;

	ifstream fin(fileName);
	string line;
	while (getline(fin, line)) {
		if (line.empty() || line[0] == '>') continue;
		refGen += line;
	}
	fin.close();

	return refGen;
}

std::map<size_t, std::pair<char, char>> FilesReader::readReferenceCsv(const string& fileName) {
	std::ifstream file(fileName);
	if (!file.is_open()) {
		std::cerr << "Failed to open the file.\n";
		throw runtime_error("Error opening file " + fileName);
	}

	std::map<size_t, std::pair<char, char>> dataMap;
	std::string line;

	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::string fStr, sStr, tStr;

		if (!std::getline(ss, fStr, ',')) continue;
		if (!std::getline(ss, sStr, ',')) continue;
		if (!std::getline(ss, tStr, ',')) continue;

		if (fStr.empty() || tStr.empty()) continue;

		//TODO remove this line
		// if (fStr[0] == 'I') continue;

		char fCol = fStr[0];
		size_t sCol = std::stoull(sStr);
		char tCol = tStr[0];

		dataMap[sCol] = std::make_pair(tCol, fCol);
	}

	return dataMap;
}

size_t FilesReader::getReferenceLength(const std::string& fileName) {
	BamReader reader;

	if (!reader.Open(fileName)) {
		std::cerr << "Error opening file " << fileName << std::endl;
		throw runtime_error("Error opening file " + fileName);
	}

	return reader.GetReferenceData()[0].RefLength;
}


AlignmentMaps FilesReader::getAlignments(const std::string& fileName, const size_t &from, const size_t &to) {
	BamReader reader;

	if (!reader.Open(fileName)) {
		std::cerr << "Error opening file " << fileName << std::endl;
		throw runtime_error("Error opening file " + fileName);
	}

	// The index file is needed in order to be able to access arbitrary reads
	if (!reader.LocateIndex() || !reader.OpenIndex(fileName + ".bai")) {
		std::cerr << "Could not open BAM index file" << fileName + ".bai" << std::endl;
		throw runtime_error("Error opening file " + fileName + ".bai");
	}

	//We are interested in reading only a certain region, not everything at once in order to not have memory exhaustion
	if (!reader.SetRegion(0, from, 0, to)) {
		std::cerr << "Error reading region" << std::endl;
	}

	//The position field in the BAM files is 0-based
	static BamAlignment alignment;
	std::map<size_t, set<pair<string, string>>> startingPos;
	std::map<size_t, pair<NucleoCounter, set<string>>> insertions;

	while (reader.GetNextAlignment(alignment)) {
		//Check whether the sequence has been aligned to the reference genome
		if (!alignment.IsMapped()) continue;

		// Clip the cigar in order to remove the clippers to be able to apply the cigar string to the expandedRead directly
		auto cigarExpanded = alignment.CigarData;
		if (cigarExpanded.begin()->Type == 'S' || cigarExpanded.begin()->Type == 'H') cigarExpanded.erase(cigarExpanded.begin());
		if ((cigarExpanded.end() - 1)->Type == 'S' || (cigarExpanded.end() - 1)->Type == 'H') cigarExpanded.erase(cigarExpanded.end() - 1);

		//AlignedBases contains the string that is clipped and has '-' inserted at the position of D in the cigar string
		string expandedRead = alignment.AlignedBases;
		//Aligned position for the string with cut out insertions (for the direct substitution and deletion analysis)
		size_t readSDStart = alignment.Position;

		//Since the insertions are cut out, if the first non-clipped part of the string requires an insertion
		//the starting index should take that into consideration
		if (cigarExpanded.begin()->Type == 'I') readSDStart += cigarExpanded.begin()->Length;
		//Base index for the insertion - should take care of the potential clipping

		size_t insertionsFound = 0;
		string noInsertionsRead;
		size_t readFromIndex = 0;
		size_t curReadIndex = 0;

		// If we analyze parts of the string that are not within the window, we can find ourselves in the situation when indices are analyzed more than once with different content
		// check whether the starting index is within the window
		// if so, iterate over the cigar string clipping it in a way that would cover only the part of the read within the window
		// store the index of the last analyzed element in the Cigar string
		// if the starting index is not within the window check the array with stored indices in order to figure out what the position that should be started with and remove the index

		for (auto iter = cigarExpanded.begin(); iter != cigarExpanded.end(); iter++) {
			//TODO ensure that only the parts of the string within the window are analyzed in order to avoid the same indices to be analyzed several times and other errors
			if (iter->Type == 'I') {
				noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
				readFromIndex = curReadIndex + iter->Length;
			}

			const bool isInsertion = iter->Type == 'I';
			for (size_t i = 0; i != iter->Length; i++) {
				const size_t insertionIndex = alignment.Position - insertionsFound;
				if (isInsertion) {
					insertions[insertionIndex + curReadIndex + i].first.increase(expandedRead[curReadIndex + i]);
					insertions[insertionIndex + curReadIndex + i].second.insert(alignment.Name);
				} else if (insertions[insertionIndex + curReadIndex + i].second.find(alignment.Name) == insertions[insertionIndex + curReadIndex + i].second.end()) {
					//Due to the difference in indexing, we may have the same sequence counted twice for the same
					insertions[insertionIndex + curReadIndex + i].first.increase('-');
				}
			}

			if (isInsertion) insertionsFound += iter->Length;
			curReadIndex += iter->Length;
		}

		//If the section before the end clipping is not an insertion one, the remaining string still needs to be read
		if (readFromIndex != curReadIndex) noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);

		// The BAM library stores the position for the aligned bases
		startingPos[readSDStart].insert(std::make_pair(noInsertionsRead, alignment.Name));
	}

	return {startingPos, insertions};
}
