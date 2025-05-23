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
	list<Read> curReads;
	list<string> curReadsNames;
	size_t linePos;
	Mutations errors;
	char curNucleo;
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

			for (linePos = 0; linePos < curRefGenLine.size(); linePos++) {
				size_t curPos = linePos + curStart;
				/*
				 * Check whether there are reads starting with pos + curPos position
				 * If so, add them to the curReads
				 * Check whether the number of current reads is smaller than 5
				 * If so, skip the analysis until
				 * Analysis:
				 *		Find the most frequent symbol for the position based on the index param of each structure
				 *		Compare it with the current one in the reference genome
				 *		If it is different, add it to the map of result denoting whether it is a deletion or a substitution
				 * Iterate over the current reads and check whether we reached the end of them
				 * If so, remove them from the list
				 * If not, increase the index param by 1
				 */

				// Add new reads that start at the current position to the list of the ones analysed
				if (startingPos.find(curPos) != startingPos.end()) {
					for (const auto& [first, second] : startingPos[curPos]) {
						curReads.emplace_front(first, curPos);
						curReadsNames.emplace_front(second);
					}
				}

				//TODO if the size of curReads is less than 5, move curPos to the next element available in startingPos or among insertion indices
				//TODO instead of constantly iterating over reads, have a map that will store pointers to the respective curReads and name of the read

				// If the number of reads for the current position is less than 5, we do not have enough data to do a meaningful evaluation
				if (curReads.size() >= MIN_READS) {
					for (auto iter = curReads.begin(); iter != curReads.end();) {
						curNucleo = (iter->sequence)[iter->index];
						nucleoCounter.increase(curNucleo);

						// If this is the end of the sequence, remove it from the list of analysed ones
						if (iter->endPos == curPos) {
							auto namesIter = curReadsNames.begin();
							advance(namesIter, distance(curReads.begin(), iter));
							curReadsNames.erase(namesIter);
							iter = curReads.erase(iter);
						} else {
							iter->index++;
							++iter;
						}
					}

					if (const char maxNucleo = nucleoCounter.findMax(curRefGenLine[linePos]); maxNucleo != curRefGenLine
						[linePos]) {
						const char actionType = maxNucleo == '-' ? 'D' : 'X';
						errors[curPos] = make_pair(maxNucleo, actionType);
					}
				} else {
					for (auto iter = curReads.begin(); iter != curReads.end();) {
						// If this is the end of the sequence, remove it from the list of analysed ones
						if (iter->endPos == curPos) {
							auto namesIter = curReadsNames.begin();
							advance(namesIter, distance(curReads.begin(), iter));
							curReadsNames.erase(namesIter);
							iter = curReads.erase(iter);
						} else {
							iter->index++;
							++iter;
						}
					}
				}

				nucleoCounter.flush();
			}
			curStart += linePos;
			linesCovered++;
		}

		errors.merge(insertions.findInsertionMutations(MIN_READS));
	}

	auto csvMap = FR::readFreeBayesVCF(referenceCsv);
	auto comp = Comparator::compareMaps(csvMap, errors);
	cout << 3;
}
