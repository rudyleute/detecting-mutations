#include <filesystem>
#include <iostream>
#include <fstream>
#include <utility>

#include "FilesManipulator.h"
#include "Comparator.h"

using FM = FilesManipulator;

#define LINES_IN_WINDOW int(1e2)
#define MIN_READS 5

using namespace std;

int main(int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();

	const string fpAlignment = FM::formFullPath(argv[1]);
	const string fpRefGen = FM::formFullPath(argv[2]);
	const string referenceCsv = FM::formFullPath(argv[3]);

	const size_t refGenLen = FM::getRefGenLength(fpAlignment);
	const string refGenName = FM::getRefGenName(fpAlignment);

	ifstream refGenFile(fpRefGen);
	if (!refGenFile.is_open()) {
		cerr << "Failed to open file: " << fpRefGen << "\n";
		return -1;
	}

	string curRefGenLine;
	size_t curStart = 0;
	Reads curReads;
	size_t linePos;
	NucleoCounter nucleoCounter;
	CompRes res;

	while (getline(refGenFile, curRefGenLine) && (curRefGenLine.empty() || curRefGenLine[0] == '>'));
	const size_t WINDOW_SIZE = curRefGenLine.size() * LINES_IN_WINDOW;

	refGenFile.seekg(0, std::ios::beg);

	auto alignments = AlignmentMaps();
	auto startingPos = alignments.startingPos;
	auto insertions = alignments.windowInsertions;
	MutationsVCF csvMap = FM::readFreeBayesVCF(referenceCsv);
	size_t windowStartInd;

	// Iterating through all the available sliding windows in order to cover the whole ref genome without memory exhaustion
	for (windowStartInd = 0; windowStartInd < refGenLen; windowStartInd += WINDOW_SIZE) {
		// Get reads within the sliding window
		alignments = FM::getAlignments(fpAlignment, windowStartInd, windowStartInd + WINDOW_SIZE, refGenName, insertions.getNextWindowInsertions());
		startingPos = alignments.startingPos;
		insertions = alignments.windowInsertions;
		Mutations errors;
		curReads.flushNonErrors();
		insertions.flushNonErrors();

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
				bool isMutation = false;
				if (csvMap.find(curPos) != csvMap.end()) {
					for (const auto& aux: csvMap.at(curPos)) {
						if (std::get<1>(aux) != 'I') {
							isMutation = true;
							break;
						}
					}
				}

				auto curErrors = curReads.iteration(curPos, linePos, MIN_READS, isMutation);
				if (!curErrors.empty()) errors.merge(curErrors);

				nucleoCounter.flush();
			}
			curStart += linePos;
			linesCovered++;
		}

		const auto insErrors = insertions.findInsertionMutations(csvMap, MIN_READS);
		for (const auto& [key, vec] : insErrors) {
			errors[key].insert(errors[key].end(), vec.begin(), vec.end());
		}
		auto nonErrors = insertions.getNonErrors();
		auto nonErrors2 = curReads.getNonErrors();
		nonErrors.insert(nonErrors2.begin(), nonErrors2.end());
		auto aux = Comparator::compareMaps(csvMap, errors, nonErrors, windowStartInd, windowStartInd + WINDOW_SIZE);
		res.merge(aux);
	}

	//We may have out of boundary indices that are not within the defined windows
	Mutations errors;
	if (auto nextWindIns = insertions.getNextWindowInsertions(); !nextWindIns.empty()) {
		errors.merge(Insertions(nextWindIns).findInsertionMutations(csvMap, MIN_READS));
		auto nonErrors = insertions.getNonErrors();
		auto nonErrors2 = curReads.getNonErrors();
		nonErrors.insert(nonErrors2.begin(), nonErrors2.end());
		auto aux = Comparator::compareMaps(csvMap, errors, nonErrors, windowStartInd, windowStartInd + WINDOW_SIZE);
		res.merge(aux);
	}

	auto end = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

	long minutes = duration.count() / 60;
	long seconds = duration.count() % 60;

	std::cout << "Execution time: " << minutes << " minutes and " << seconds << " seconds" << std::endl;

	FilesManipulator::saveToCsv(FM::getRefGenName(fpAlignment) + "new", res);
}
