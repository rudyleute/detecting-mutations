#ifndef COMPARATOR_H
#define COMPARATOR_H

#include <map>
#include <string>
#include <vector>

#include "FilesReader.h"

using namespace std;

struct compRes;
struct Read;

class Comparator {
public:
	static compRes compareMaps(
		const Mutations& map1,
		const Mutations& map2
	);
};

struct compRes {
	std::vector<std::tuple<size_t, char, char>> difference1;
	std::vector<std::tuple<size_t, char, char>> difference2;
	std::vector<std::tuple<size_t, char, char, char, char>> errors;

	explicit compRes(
		std::vector<std::tuple<size_t, char, char>> difference1,
		std::vector<std::tuple<size_t, char, char>> difference2,
		vector<std::tuple<size_t, char, char, char, char>> errors
	): difference1(std::move(difference1)), difference2(std::move(difference2)), errors(std::move(errors)) {}
};

struct Read {
	size_t index;
	const size_t endPos;
	const string sequence;

	explicit Read(const string& read, const size_t& startingPos): index(0),
	                                                              endPos(startingPos + read.size() - 1),
	                                                              sequence(read) {}
};

struct Reads {
private:
	list<Read> reads;
	list<string> names;
	string curRefGenLine;

public:
	Reads() = default;

	void setRefGenLine(const string& refGenLine) {
		this->curRefGenLine = refGenLine;
	}

	void addReads(const Alignments& startingPos, const size_t &curPos) {
		for (const auto& [first, second] : startingPos.at(curPos)) {
			reads.emplace_front(first, curPos);
			names.emplace_front(second);
		}
	}

	Mutations iteration(const size_t& curPos, const size_t& linePos, const size_t &minReads) {
		char curNucleo;
		NucleoCounter nucleoCounter;
		Mutations errors;

		const bool isRelevant = reads.size() >= minReads;
		for (auto iter = reads.begin(); iter != reads.end();) {
			if (isRelevant) {
				curNucleo = (iter->sequence)[iter->index];
				nucleoCounter.increase(curNucleo);
			}

			// If this is the end of the sequence, remove it from the list of analysed ones
			if (iter->endPos == curPos) {
				auto namesIter = names.begin();
				advance(namesIter, distance(reads.begin(), iter));
				names.erase(namesIter);
				iter = reads.erase(iter);
			} else {
				iter->index++;
				++iter;
			}
		}

		if (isRelevant) {
			if (const char maxNucleo = nucleoCounter.findMax(curRefGenLine[linePos]); maxNucleo != curRefGenLine[linePos]) {
				const char actionType = maxNucleo == '-' ? 'D' : 'X';
				errors[curPos] = make_pair(maxNucleo, actionType);
			}
		}

		return errors;
	}
};


#endif
