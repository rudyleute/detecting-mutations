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
	static compRes compareMaps(const Mutations& map1,
	                    const Mutations& map2);
};

struct compRes {
	std::vector<std::tuple<size_t, char, char>> difference1;
	std::vector<std::tuple<size_t, char, char>> difference2;
	std::vector<std::tuple<size_t, char, char, char, char>> errors;

	explicit compRes(std::vector<std::tuple<size_t, char, char>> difference1,
	                 std::vector<std::tuple<size_t, char, char>> difference2,
	                 vector<std::tuple<size_t, char, char, char, char>> errors):
		difference1(std::move(difference1)), difference2(std::move(difference2)), errors(std::move(errors)) {}
};

struct Read {
	size_t index;
	const size_t endPos;
	const string sequence;

	explicit Read(const string& read, const size_t& startingPos): index(0), endPos(startingPos + read.size() - 1),
																  sequence(read) {}
};


#endif
