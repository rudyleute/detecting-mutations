#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

struct CompRes;
struct NucleoCounter;
struct Read;
struct AlignmentMaps;
struct Insertions;

using namespace std;
using MutationErrors = std::vector<std::tuple<size_t, char, char, NucleoCounter>>;
using InBothEr = std::vector<std::tuple<size_t, char, char, char, char, NucleoCounter>>;
using InsertionMap = std::map<size_t, std::pair<NucleoCounter, std::set<std::string>>>;
using CigarString = list<pair<char, size_t>>;
using Mutations = std::map<size_t, std::tuple<char, char, NucleoCounter>>;
using MutationsVCF = std::map<size_t, std::tuple<char, char>>;
using Alignments = map<size_t, set<pair<string, string>>>;

inline std::map<char, size_t> nucleoMapping = {
	{'-', 0},
	{'A', 1},
	{'G', 2},
	{'C', 3},
	{'T', 4},
};

inline std::map<size_t, char> rNucleoMapping = [] {
	std::map<size_t, char> reverse;
	for (const auto& [key, value] : nucleoMapping) {
		reverse[value] = key;
	}
	return reverse;
}();

struct NucleoCounter {
private:
	std::vector<size_t> counters;

public:
	NucleoCounter(): counters(nucleoMapping.size()) {}

	void increase(const char& nucleo) {
		counters[nucleoMapping[nucleo]]++;
	}

	void setCounter(const char& nucleo, const size_t& value) {
		counters[nucleoMapping[nucleo]] = value;
	}

	char findMax(const char& base) const {
		size_t curMax = 0;
		set<char> curMaxNucleos;

		for (size_t i = 0; i != counters.size(); i++) {
			if (counters[i] >= curMax) {
				if (counters[i] > curMax) {
					curMax = counters[i];
					curMaxNucleos.clear();
				}

				curMaxNucleos.insert(rNucleoMapping[i]);
			}
		}

		const double ratio = static_cast<double>(curMax) / size();
		//TODO probably it is worth considering adding an epsilon here?
		if (ratio >= 0.5) {
			//Return the first alphabetically sorted non-base value if it is seen in at least 50% of cases
			if (ratio > 0.5 || curMaxNucleos.size() == 1) return *curMaxNucleos.begin();
			for (const auto& elem : curMaxNucleos) if (elem != base) return elem;
		}

		return base;
	}

	void flush() {
		std::fill(counters.begin(), counters.end(), 0);
	}

	void merge(const NucleoCounter& other) {
		const auto result = other.getCounters();
		for (size_t i = 0; i != counters.size(); i++) {
			counters[i] += result[i];
		}
	}

	vector<size_t> getCounters() const {
		return counters;
	}

	size_t size() const {
		size_t counter = 0;
		for (const size_t i : counters) counter += i;

		return counter;
	}
};

struct Insertions {
private:
	InsertionMap insertions;
	InsertionMap nextWindowInsertions;
	set<size_t> insertionIndices;

	InsertionMap* curMap;
	string expandedRead;
	string name;

	std::map<size_t, NucleoCounter> nonErrors;

	void addValues(const size_t& refGenIndex, const size_t& curReadIndex, const size_t& start, const size_t& end, const bool isInsertion) {
		for (size_t i = start; i != end; i++) {
			if (isInsertion) {
				(*curMap)[refGenIndex + i].first.increase(expandedRead[curReadIndex + i]);
				(*curMap)[refGenIndex + i].second.insert(name);
				insertionIndices.insert(refGenIndex + i);
			} else if ((*curMap)[refGenIndex + i].second.find(name) == (*curMap)[refGenIndex + i].second.end()) {
				(*curMap)[refGenIndex + i].first.increase('-');
			}
		}
	}

public:
	Insertions() = default;

	Insertions(
		const InsertionMap& insertions,
		const InsertionMap& nextWindowInsertions = InsertionMap(),
		const set<size_t>& insertionIndices = set<size_t>()
	): insertions(insertions), nextWindowInsertions(nextWindowInsertions), insertionIndices(insertionIndices) {
		curMap = &this->insertions;
	}

	void setRead(const string& expandedRead, const string& name) {
		this->expandedRead = move(expandedRead);
		this->name = move(name);
	}
	InsertionMap getNextWindowInsertions() const {
		return nextWindowInsertions;
	}

	void addInsertion(
		const size_t& refGenIndex,
		const size_t& curReadIndex,
		const size_t& end,
		const size_t& left,
		const bool isInsertion = false
	) {
		if (left != end) {
			this->addValues(refGenIndex, curReadIndex, 0, left, isInsertion);
			curMap = &this->nextWindowInsertions;
			this->addValues(refGenIndex, curReadIndex, left, end, isInsertion);
			curMap = &this->insertions;
		} else this->addValues(refGenIndex, curReadIndex, 0, end, isInsertion);
	}

	Mutations findInsertionMutations(const MutationsVCF &mutationsVCF, const size_t &minReads) {
		Mutations errors;

		for (const size_t &index: insertionIndices) {
			if (insertions[index].first.size() >= minReads)
				if (const char maxNucleo = insertions[index].first.findMax('-'); maxNucleo != '-') errors[index] = make_tuple(
					maxNucleo, 'I', insertions[index].first);
				else if (mutationsVCF.find(index) != mutationsVCF.end() && std::get<1>(mutationsVCF.at(index)) == 'I') {
					nonErrors[index] = insertions[index].first;
				}
		}

		return errors;
	}

	std::map<size_t, NucleoCounter> getNonErrors() {
		return nonErrors;
	}

	void flushNonErrors() {
		nonErrors.clear();
	}
};

struct AlignmentMaps {
	Alignments startingPos;
	Insertions windowInsertions;

	AlignmentMaps(
		const Alignments& alignments,
		const Insertions& windowInsertions
	): startingPos(alignments),
	   windowInsertions(windowInsertions) {}

	AlignmentMaps() = default;
};

struct CompRes {
	MutationErrors diffInVCF;
	MutationErrors diffInCust;
	InBothEr errors;
	std::map<size_t, NucleoCounter> nonErrors;


	CompRes() = default;
	explicit CompRes(
		MutationErrors diffInVCF,
		MutationErrors diffInCust,
		InBothEr errors
	): diffInVCF(std::move(diffInVCF)), diffInCust(std::move(diffInCust)), errors(std::move(errors)) {}

	void merge(CompRes &res) {
		diffInVCF.insert(diffInVCF.end(), res.diffInVCF.begin(), res.diffInVCF.end());
		diffInCust.insert(diffInCust.end(), res.diffInCust.begin(), res.diffInCust.end());
		this->errors.insert(this->errors.end(), res.errors.begin(), res.errors.end());
	}
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
	std::map<size_t, NucleoCounter> nonErrors;

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

	Mutations iteration(const size_t& curPos, const size_t& linePos, const size_t &minReads, const bool isReported) {
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
				errors[curPos] = make_tuple(maxNucleo, actionType, nucleoCounter);
			} else if (isReported) {
				nonErrors[curPos] = nucleoCounter;
			}
		}

		return errors;
	}

	std::map<size_t, NucleoCounter> getNonErrors() {
		return nonErrors;
	}

	void flushNonErrors() {
		nonErrors.clear();
	}
};

#endif //STRUCTURES_H
