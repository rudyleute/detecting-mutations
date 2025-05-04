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
		if (fStr[0] == 'I') continue;

		char fCol = fStr[0];
		size_t sCol = std::stoull(sStr);
		char tCol = tStr[0];

		dataMap[sCol] = std::make_pair(tCol, fCol);
	}

	return dataMap;
}

AlignmentMaps FilesReader::getAlignments(const std::string& fileName) {
	BamReader reader;

	if (!reader.Open(fileName)) {
		std::cerr << "Error opening file " << fileName << std::endl;
		throw runtime_error("Error opening file " + fileName);
	}

	//The position field in the BAM files is 0-based
	static BamAlignment alignment;
	std::map<size_t, set<string>> startingPos;
	std::map<size_t, std::map<char, int>> insertions;

	while (reader.GetNextAlignment(alignment)) {
		if (!alignment.IsMapped()) continue;
		string expandedRead = alignment.AlignedBases;

		//The CigarData is constructed for the non-clipped string and we are analyzing the clipped one
		size_t clippingOffset = 0;
		if ((alignment.CigarData.begin())->Type == 'S') clippingOffset = (alignment.CigarData.begin())->Length;
		//If the first non-clipping action is insertion, then it will be cut out and the starting point of the string changes
		//The addition of the insertion symbols will be performed within the next cycle
		size_t offset = 0;
		if ((alignment.CigarData.begin() + 1)->Type == 'I') offset = (alignment.CigarData.begin() + 1)->Length;

		string noInsertionsRead;
		size_t readFromIndex = 0;
		size_t curIndex = -clippingOffset;
		for (auto iter = alignment.CigarData.begin(); iter != alignment.CigarData.end(); iter++) {
			if (iter->Type == 'I') {
				noInsertionsRead += expandedRead.substr(readFromIndex, curIndex - readFromIndex);
				readFromIndex = curIndex + iter->Length;

				for (size_t i = 0; i != iter->Length; i++) insertions[curIndex + i][expandedRead[curIndex + i]]++;
			}

			curIndex += iter->Length;
		}
		if ((alignment.CigarData.end() - 1)->Type == 'S') curIndex -= (alignment.CigarData.end() - 1)->Length;
		if (readFromIndex != curIndex) noInsertionsRead += expandedRead.substr(readFromIndex, curIndex - readFromIndex);

		// The BAM library stores the position for the aligned bases
		startingPos[alignment.Position + offset].insert(noInsertionsRead);

		// //#TODO ensure that this is a correct way to calculate the insertion indices
		// size_t curPos = alignment.Position;
		// for (auto iterator = alignment.CigarData.begin(); iterator != alignment.CigarData.end() - 1; ++iterator) {
		// 	if (iterator->Type == 'I') {
		// 		for (size_t insertionInd = curPos; insertionInd < curPos + iterator->Length; insertionInd++) {
		// 			insertionMap[insertionInd]++;
		// 		}
		// 	}
		// 	curPos += iterator->Length;
		// }
	}

	return {startingPos, insertions};
}
