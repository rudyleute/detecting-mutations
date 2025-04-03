#include <iostream>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include "SamReader.h"
#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval.hpp>
#include <map>
#include <utility>
#include <vector>
#include <string>

using namespace BamTools;

char SamReader::complementSymbol(const char& n) {
	switch (n) {
		case 'A': return 'T';
		case 'T': return 'A';
		case 'G': return 'C';
		case 'C': return 'G';
		case '-': return '-';
		default: return 'N';
	}
}

void SamReader::reverseComplement(std::string& read) {
	reverse(read.begin(), read.end());
	transform(read.begin(), read.end(), read.begin(), complementSymbol);
}

AlignmentMaps SamReader::getAlignments(const std::string& fileName) {
	BamReader reader;

	if (!reader.Open(fileName)) {
		std::cerr << "Error opening file " << fileName << std::endl;
		throw runtime_error("Error opening file " + fileName);
	}

	//The position field in the BAM and SAM files is 1-based
	static BamAlignment alignment;
	interval_map<size_t, string> intervalAlignments;
	std::map<size_t, size_t> insertionMap;

	while (reader.GetNextAlignment(alignment)) {
		if (!alignment.IsMapped()) continue;

		string aln = alignment.AlignedBases;
		if (alignment.IsReverseStrand()) SamReader::reverseComplement(aln);

		size_t endPos = alignment.Position + aln.length();
		intervalAlignments.insert(make_pair(interval<size_t>::right_open(alignment.Position, endPos), aln));

		//#TODO ensure that this is a correct way to calculate the insertion indices
		size_t curPos = alignment.Position;
		for (auto iterator = alignment.CigarData.begin(); iterator != alignment.CigarData.end() - 1; ++iterator) {
			if (iterator->Type == 'I') {
				for (size_t insertionInd = curPos; insertionInd < curPos + iterator->Length; insertionInd++) {
					insertionMap[insertionInd]++;
				}
			}
			curPos += iterator->Length;
		}
	}

	return AlignmentMaps(intervalAlignments, insertionMap);
}
