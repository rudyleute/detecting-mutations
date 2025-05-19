#include <filesystem>
#include <iostream>
#include <fstream>
#include <utility>
#include <list>

#include "FilesReader.h"

#define FASTA_LINE_LEN 80
#define LINES_IN_WINDOW int(1e2)
#define WINDOW_SIZE (FASTA_LINE_LEN * LINES_IN_WINDOW)
#define MIN_READS 5

using namespace std;

struct Read {
	size_t index;
	const size_t endPos;
	const string sequence;

	explicit Read(const string& read, const size_t& startingPos): index(0), endPos(startingPos + read.size() - 1),
	                                                              sequence(read) {}
};

struct compRes {
	std::vector<std::tuple<size_t, char, char>> difference1;
	std::vector<std::tuple<size_t, char, char>> difference2;
	std::vector<std::tuple<size_t, char, char, char, char>> errors;

	explicit compRes(std::vector<std::tuple<size_t, char, char>> difference1,
	                 std::vector<std::tuple<size_t, char, char>> difference2,
	                 std::vector<std::tuple<size_t, char, char, char, char>> errors):
		difference1(std::move(difference1)), difference2(std::move(difference2)), errors(std::move(errors)) {}
};

compRes compareMaps(const std::map<size_t, std::pair<char, char>>& map1,
                    const std::map<size_t, std::pair<char, char>>& map2) {
	std::vector<std::tuple<size_t, char, char>> difference1;
	std::vector<std::tuple<size_t, char, char>> difference2;
	std::vector<std::tuple<size_t, char, char, char, char>> errors;

	for (const auto& [key, val1] : map1) {
		auto it = map2.find(key);
		if (it == map2.end()) difference1.emplace_back(key, val1.first, val1.second);
		else if (it->second != val1)
			errors.emplace_back(key, val1.first, val1.second, it->second.first,
			                    it->second.second);
	}

	for (const auto& [key, val2] : map2) {
		if (map1.find(key) == map1.end()) difference2.emplace_back(key, val2.first, val2.second);
	}

	return compRes(difference1, difference2, errors);
}

size_t findDifferenceSize(const list<string>& first, const set<string>& second) {
	set<string> result;

	for (const auto& elem : first)
		if (second.find(elem) != second.end()) result.insert(elem);

	return result.size();
}

int main(int argc, char* argv[]) {
	const std::string fpAlignment = std::filesystem::current_path().parent_path().string() + '/' + argv[1];
	const std::string fpRefGen = std::filesystem::current_path().parent_path().string() + '/' + argv[2];
	const std::string referenceCsv = std::filesystem::current_path().parent_path().string() + '/' + argv[3];

	auto refGenLen = FilesReader::getReferenceLength(fpAlignment);

	ifstream refGenFile(fpRefGen);
	if (!refGenFile.is_open()) {
		std::cerr << "Failed to open file: " << fpRefGen << "\n";
		return -1;
	}

	string line;
	size_t curStart = 0;
	list<Read> curReads;
	list<string> curReadsNames;
	size_t linePos;
	std::map<size_t, std::pair<char, char>> errors;
	char curNucleo;
	NucleoCounter nucleoCounter;

	// Iterating through all the available sliding windows in order to cover the whole ref genome without memory exhaustion
	for (size_t windowStartInd = 0; windowStartInd < refGenLen; windowStartInd += WINDOW_SIZE) {
		// Get reads within the sliding window
		auto alignments = FilesReader::getAlignments(fpAlignment, windowStartInd, windowStartInd + WINDOW_SIZE);
		auto startingPos = alignments.startingPos;
		auto insertions = alignments.insertionCount;

		size_t linesCovered = 0;
		while (linesCovered < LINES_IN_WINDOW && getline(refGenFile, line)) {
			if (line.empty() || line[0] == '>') continue;

			for (linePos = 0; linePos < line.size(); linePos++) {
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
				if (startingPos.find(curPos) != startingPos.end())
					for (const auto& [first, second] : startingPos[curPos]) {
						curReads.emplace_front(first, curPos);
						curReadsNames.emplace_front(second);
					}

				//TODO if the size of curReads is less than 5, move curPos to the next element available in startingPos or among insertion indices

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
						}
						else {
							iter->index++;
							++iter;
						}
					}

					if (const char maxNucleo = nucleoCounter.findMax(line[linePos]); maxNucleo != line[linePos]) {
						const char actionType = maxNucleo == '-' ? 'D' : 'X';
						errors[curPos] = std::make_pair(maxNucleo, actionType);
					}
				} else {
					for (auto iter = curReads.begin(); iter != curReads.end();) {
						// If this is the end of the sequence, remove it from the list of analysed ones
						if (iter->endPos == curPos) {
							auto namesIter = curReadsNames.begin();
							advance(namesIter, distance(curReads.begin(), iter));
							curReadsNames.erase(namesIter);
							iter = curReads.erase(iter);
						}
						else {
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

		//TODO ensure that only indices within the current window range are analyzed. The rest must be propagated to the next iteration
		for (const auto& [first, second] : insertions) {
			if (first == 846) {
				cout << 2;
			}
			if (second.first.size() >= MIN_READS)
				if (const char maxNucleo = second.first.findMax('-', true); maxNucleo != '-')
					errors[first] = std::make_pair(maxNucleo, 'I');
		}
	}

	auto csvMap = FilesReader::readReferenceCsv(referenceCsv);
	auto comp = compareMaps(csvMap, errors);
	cout << 2 << std::endl;
}
