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

		//AlignedBases contains the string that is clipped and has '-' inserted at the position of D in the cigar string
		string expandedRead = alignment.AlignedBases;
		//Aligned position for the string with cut out insertions (for the direct substitution and deletion analysis)
		size_t readSDStart = alignment.Position;

		//Current index in the read - we are traversing the read and use it. Can be negative at first if there is a clipping
		size_t curReadIndex = 0;
		if ((alignment.CigarData.begin())->Type == 'S' || (alignment.CigarData.begin())->Type == 'H') curReadIndex -= (alignment.CigarData.begin())->Length;

		//Since the insertions are cut out, if the first non-clipped part of the string requires an insertion
		//the starting index should take that into consideration
		if ((alignment.CigarData.begin() + 1)->Type == 'I') readSDStart += (alignment.CigarData.begin() + 1)->Length;
		//Base index for the insertion - should take care of the potential clipping
		size_t insertionBaseIndex = alignment.Position;

		string noInsertionsRead;
		size_t readFromIndex = 0;
		for (auto iter = alignment.CigarData.begin(); iter != alignment.CigarData.end(); iter++) {
			//TODO ensure that only the parts of the string within the window are analyzed in order to avoid the same indices to be analyzed several times and other errors
			if (iter->Type == 'I') {
				noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
				readFromIndex = curReadIndex + iter->Length;
			}

			if (iter->Type != 'S') {
				const bool isInsertion = iter->Type == 'I';
				for (size_t i = 0; i != iter->Length; i++) {
					if (isInsertion) {
						insertions[insertionBaseIndex + curReadIndex + i].first.increase(expandedRead[curReadIndex + i]);
						insertions[insertionBaseIndex + curReadIndex + i].second.insert(alignment.Name);
					} else if (insertions[insertionBaseIndex + curReadIndex + i].second.find(alignment.Name) == insertions[insertionBaseIndex + curReadIndex + i].second.end()) {
						//Due to the difference in indexing, we may have the same sequence counted twice for the same
						insertions[insertionBaseIndex + curReadIndex + i].first.increase('-');
					}
				}

				if (isInsertion) insertionBaseIndex -= iter->Length;
			}

			curReadIndex += iter->Length;
		}

		//If clipping is needed at the end of the string as well, then in the cycle above the index is off
		if ((alignment.CigarData.end() - 1)->Type == 'S' || (alignment.CigarData.end() - 1)->Type == 'H') curReadIndex -= (alignment.CigarData.end() - 1)->Length;
		//If the section before the end clipping is not an insertion one, the remaining string still needs to be read
		if (readFromIndex != curReadIndex) noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);

		// The BAM library stores the position for the aligned bases
		startingPos[readSDStart].insert(std::make_pair(noInsertionsRead, alignment.Name));
	}

	return {startingPos, insertions};
}
