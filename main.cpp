#include <filesystem>
#include <iostream>
#include <fstream>
#include <utility>
#include <list>

#include "FilesReader.h"
#include "Comparator.h"

using FR = FilesReader;

#define FASTA_LINE_LEN 80
#define LINES_IN_WINDOW int(1e2)
#define WINDOW_SIZE (FASTA_LINE_LEN * LINES_IN_WINDOW)
#define MIN_READS 5

using namespace std;

int main(int argc, char* argv[]) {
	const string fpAlignment = FR::formFullPath(argv[1]);
	const string fpRefGen = FR::formFullPath(argv[2]);
	const string referenceCsv = FR::formFullPath(argv[3]);

	auto refGenLen = FR::getReferenceLength(fpAlignment);

	ifstream refGenFile(fpRefGen);
	if (!refGenFile.is_open()) {
		cerr << "Failed to open file: " << fpRefGen << "\n";
		return -1;
	}

	string curRefGenLine;
	size_t curStart = 0;
	Reads curReads;
	size_t linePos;
	Mutations errors;
	NucleoCounter nucleoCounter;

	auto alignments = AlignmentMaps();
	auto startingPos = alignments.startingPos;
	auto insertions = alignments.windowInsertions;

	// Iterating through all the available sliding windows in order to cover the whole ref genome without memory exhaustion
	for (size_t windowStartInd = 0; windowStartInd < refGenLen; windowStartInd += WINDOW_SIZE) {
		// Get reads within the sliding window
		alignments = FR::getAlignments(fpAlignment, windowStartInd, windowStartInd + WINDOW_SIZE, insertions.getNextWindowInsertions());
		startingPos = alignments.startingPos;
		insertions = alignments.windowInsertions;

		size_t linesCovered = 0;
		while (linesCovered < LINES_IN_WINDOW && getline(refGenFile, curRefGenLine)) {
			if (curRefGenLine.empty() || curRefGenLine[0] == '>') continue;

			curReads.setRefGenLine(curRefGenLine);
			for (linePos = 0; linePos < curRefGenLine.size(); linePos++) {
				size_t curPos = linePos + curStart;

				// Add new reads that start at the current position to the list of the ones analysed
				if (startingPos.find(curPos) != startingPos.end()) curReads.addReads(startingPos, curPos);

				//TODO if the size of curReads is less than 5, move curPos to the next element available in startingPos or among insertion indices
				//TODO instead of constantly iterating over reads, have a map that will store pointers to the respective curReads and name of the read

				// If the number of reads for the current position is less than 5, we do not have enough data to do a meaningful evaluation
				auto curErrors = curReads.iteration(curPos, linePos, MIN_READS);
				if (!curErrors.empty()) errors.merge(curErrors);

				nucleoCounter.flush();
			}
			curStart += linePos;
			linesCovered++;
		}

		errors.merge(insertions.findInsertionMutations(MIN_READS));
	}

	//We may have out of boundary indices that are not within the defined windows
	if (auto nextWindIns = insertions.getNextWindowInsertions(); !nextWindIns.empty()) {
		errors.merge(Insertions(nextWindIns).findInsertionMutations(MIN_READS));
	}

	auto csvMap = FR::readFreeBayesVCF(referenceCsv);
	auto comp = Comparator::compareMaps(csvMap, errors);
	cout << 3;
}
